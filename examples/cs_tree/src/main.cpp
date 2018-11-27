//#define VERBOSE

#include <iostream>
#include <fstream>
#include <cstring>
#include <list>
#include <cassert>
#include <PiecewisePolynomialFunction.h>
#include <RoundingClasses.h>
#include <ContourSpectrum.h>
#include <sys/types.h>

extern "C" 
{
#include <tourtre.h>
}

#include "Data.h"
#include "Mesh.h"

#include <unistd.h>

// Contour tree calculation helpers
double value (size_t v, void * d)
{
    Mesh * mesh = reinterpret_cast<Mesh*>(d);
    return mesh->data[v];
}


size_t neighbors ( size_t v, size_t * nbrs, void * d )
{
    Mesh *mesh = static_cast<Mesh*>(d);
    static std::vector<size_t> nbrsBuf;

    nbrsBuf.clear();
    mesh->getNeighbors(v, nbrsBuf);

    for (uint i = 0; i < nbrsBuf.size(); i++) {
        nbrs[i] = nbrsBuf[i]; 
    }
    return nbrsBuf.size();
}

// Data structures
typedef int IdT;
#define ID_TO_VOIDPTR(id) reinterpret_cast<void*>(id)
#define VOIDPTR_TO_ID(ptr) static_cast<IdT>(reinterpret_cast<uint64_t>(ptr))

IdT currentId;

struct BranchData
{
    int depth;
    IdT id;
    PiecewisePolynomialFunction<double> areaPPF;
    PiecewisePolynomialFunction<double> volumePPF;
};

#define ITERATE_OVER_CHILDREN(it, b) \
    for (ctBranch* it=(b)->children.head; it != NULL; it = it->nextChild)

// Add id to branch. Abuse data pointer to store id
IdT addId(ctBranch* b, IdT currId=0)
{
    b->data = ID_TO_VOIDPTR(currId++);

    ITERATE_OVER_CHILDREN(c, b)
    {
        currId = addId(c, currId);
    }

    return currId;
}

// Convert a branch map that has ids instead of pointers to branches
IdT* createIdBasedBranchMap(ctBranch **map, size_t mapSize)
{
    IdT* idBasedMap = new IdT[mapSize];

    // Iterate over pointer based map and read id for each branch
    for (int idx=0; idx<mapSize; ++idx)
        idBasedMap[idx] = VOIDPTR_TO_ID(map[idx]->data);

    return idBasedMap;
}

// Add id to translation map for all children. This function is used
// for pruned away portions of the tree.
void addEntriesToOldIdToNewIdMap(ctBranch *b, IdT* oldIdToNewIdMap, IdT& newId)
{
    oldIdToNewIdMap[VOIDPTR_TO_ID(b->data)] = newId;
    ITERATE_OVER_CHILDREN(c, b)
    {
        addEntriesToOldIdToNewIdMap(c, oldIdToNewIdMap, newId);
    }
}

// Copy tree to a version that removes symbolic perturbation. Also add branch data and 
// create a map that translates between ids from the unsimplified to the simplified
// branch decomposition.
ctBranch* ctBranch_copyRemovePerturbationAndAddBranchData(ctBranch *b, Data& d, ctContext* ctx, IdT* oldIdToNewIdMap, IdT& currId, int depth = 0)
{
    // Create copy of current branch
    ctBranch *res = ctBranch_new(b->extremum, b->saddle, ctx);

    // Add branch data
    BranchData *bd = new BranchData;
    res->data = bd;
    bd->depth = depth;
    bd->id = currId++;

    // Add translation map entry
    oldIdToNewIdMap[VOIDPTR_TO_ID(b->data)] = bd->id;

    ITERATE_OVER_CHILDREN(c, b)
    {
        if (fabs(double(d[c->extremum]) - double(d[c->saddle])) > 0)
        {
            ctBranch *newC = ctBranch_copyRemovePerturbationAndAddBranchData(c, d, ctx, oldIdToNewIdMap, currId, depth + 1);
            ctBranchList_add( &(res->children), newC , ctx);
            newC->parent = res;
        }
        else 
        {
            addEntriesToOldIdToNewIdMap(c, oldIdToNewIdMap, bd->id);
        }
    }
    return res;
}

void applyIdTranslationToBranchMap(IdT* map, size_t mapSize, IdT* oldIdToNewIdMap)
{
    for (int idx=0; idx<mapSize; ++idx)
        map[idx] = oldIdToNewIdMap[map[idx]];
}

// Sanity check for simplified tree
void checkTree(ctBranch *b, Data& d)
{
    // All zero persistence branches need to be simplified away
    assert(d[b->extremum] != d[b->saddle]);
    // All branches except root need to have non-zero parent pointer
    assert(static_cast<BranchData*>(b->data)->depth == 0 || b->parent);

    ITERATE_OVER_CHILDREN(c, b)
        checkTree(c, d);
}

void writeIdIntoMap(ctBranch** map, ctBranch *b)
{
    map[static_cast<BranchData*>(b->data)->id] = b;
    ITERATE_OVER_CHILDREN(c, b)
        writeIdIntoMap(map, c);
}

ctBranch** createIdToBranchPtrMap(ctBranch*b, size_t mapSize)
{
    ctBranch** newMap = new ctBranch*[mapSize];
    for (int idx = 0; idx < mapSize; ++idx)
        newMap[idx]=0;
    writeIdIntoMap(newMap, b);
    for (int idx = 0; idx<mapSize; ++idx)
        assert(newMap[idx]);
    return newMap;
}

double computeVolumeFromTree(ctBranch *b, Data& d, double val)
{
    double eV = d[b->extremum];
    double sV = d[b->saddle];
    double mV = std::max(eV, sV);
    PiecewisePolynomialFunction<double>& volPPF = static_cast<BranchData*>(b->data)->volumePPF;
    double vol = 0;
    std::list<double> saddleVals;
    ITERATE_OVER_CHILDREN(c, b)
    {
        if (d[c->saddle] <= val) saddleVals.push_back(d[c->saddle]);
        vol += computeVolumeFromTree(c, d, val);
    }
    if (!volPPF.empty())
    {
        saddleVals.sort();
        saddleVals.unique();
        for (std::list<double>::iterator sit = saddleVals.begin(); sit!=saddleVals.end(); ++ sit)
           vol += volPPF(*sit-0.000000001);
        vol += volPPF(std::min(val, mV));
    }
    return vol;
}

double computeAreaFromTree(ctBranch *b, Data& d, double val)
{
    PiecewisePolynomialFunction<double>& areaPPF = static_cast<BranchData*>(b->data)->areaPPF;
    double area =0;
    ITERATE_OVER_CHILDREN(c, b)
        area += computeAreaFromTree(c, d, val);
    if (!areaPPF.empty())
        area += areaPPF(val);
    return area;
}

void savePPF(PiecewisePolynomialFunction<double>& ppf, FILE *f)
{
    std::vector<char> stream;
    std::back_insert_iterator<std::vector<char> > ins(std::back_inserter(stream));
    ppf.toByteStream(ins);
    size_t ppfSize = stream.size();
    fwrite(&ppfSize, sizeof(size_t), 1, f);
    if (ppfSize > 0)
        fwrite(&stream[0], 1, stream.size(), f);
    else
        std::cout << "Warning: Writing empty piecewise linear function!" << std::endl;
}

void saveCS(void *d, FILE *f)
{
    BranchData *bd = static_cast<BranchData*>(d);
    if (bd)
    {
        size_t numSig = 2;
        fwrite(&numSig, sizeof(size_t), 1, f);
        savePPF(bd->areaPPF, f);
        savePPF(bd->volumePPF, f);
    }
    else
    {
        size_t numSig = 0;
        fwrite(&numSig, sizeof(size_t), 1, f);
    }
}

void outputTree(std::ostream & out, ctBranch * b, Data &d)
{
    if (fabs(double(d[b->extremum]) - double(d[b->saddle])) > 0)
    {
	out << "(" << int(d[b->extremum]) << " [" << b->extremum << "] " << int(d[b->saddle]) << " [" << b->saddle << "]";
        if (b->data)
        {
            out << " depth = " << static_cast<BranchData*>(b->data)->depth << " id = " << static_cast<BranchData*>(b->data)->id;
            //out << " area = ";
            //static_cast<BranchData*>(b->data)->areaPPF.print(out);
            //char filename[1024];
            //sprintf(filename,"volume_%d_%d_%d_%d.dat", int(d[b->extremum]), int(b->extremum), int(d[b->saddle]), int(b->saddle));
            //std::cout << "Saving curve " << filename << std::endl;
            //static_cast<BranchData*>(b->data)->volumePPF.saveCurve(filename,200);
        }
	
        ITERATE_OVER_CHILDREN(c, b)
        {
		out << " ";
		outputTree( out, c, d );
	}
	
	out << ")";
    }
} 

// Find the path between two points in the contour tree
void findPath(ctBranch* b1, ctBranch* b2, std::list<ctBranch*> &fromB1, std::list<ctBranch*> &fromB2, ctBranch *&commonParentBranch)
{
    fromB1.clear();
    fromB2.clear();

    assert(b1->data && b2->data);
    int d1 = static_cast<BranchData*>(b1->data)->depth;
    int d2 = static_cast<BranchData*>(b2->data)->depth;

    if (d1 > d2)
    {
        while (d1 > d2)
        {
            //assert(b1->parent);
            fromB1.push_back(b1);
            b1 = b1->parent;
            --d1;
            //assert(d1 == static_cast<BranchData*>(b1->data)->depth);
        }
    }
    else if (d2 > d1)
    {
        while (d2 > d1)
        {
            //assert(b2->parent);
            fromB2.push_front(b2);
            b2 = b2->parent;
            --d2;
            //assert(d2 == static_cast<BranchData*>(b2->data)->depth);
        }
    }

    //assert (d1 == d2);

    while (b1 != b2)
    {
        //assert(b1->parent);
        //assert(b2->parent);
        fromB1.push_back(b1);
        b1=b1->parent;
        fromB2.push_front(b2);
        b2=b2->parent;
        //assert(static_cast<BranchData*>(b1->data)->depth == static_cast<BranchData*>(b2->data)->depth);
        //assert(b1);
        //assert(b2);
    }
    commonParentBranch = b1;
}

void calcContourSpectrum(Data &d, ctBranch** bm)
{
    RoundToNearestMultiple roundF(0.1);
    ContourSpectrum<RoundToNearestMultiple> cs(true, true, false, roundF);

    for (unsigned int i=0; i<d.size[0]-1; ++i)
    {
        for (unsigned int j=0; j<d.size[1]-1; ++j)
            for (unsigned int k=0; k<d.size[2]-1; ++k)
            {
#ifdef VERBOSE
                std::cout << "i="<<i<<" j="<<j<<" k="<<k << std::endl;
#endif
                unsigned int pos[8][3] = {
                    { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 },
                    { 0, 0, 1 }, { 1, 0, 1 }, { 1, 1, 1 }, { 0, 1, 1  } };

                unsigned int ptId[8];
                for (int vtx=0; vtx<8; ++vtx)
                    ptId[vtx] = d.convertIndex(i+pos[vtx][0], j+pos[vtx][1], k+pos[vtx][2]);

                unsigned int tets[2][5][4] = {
                    { {0, 2, 5, 1}, {0, 7, 2, 3}, {0, 5, 7, 4}, {2, 7, 5, 6}, {0, 7, 5, 2} },
                    { {1, 4, 3, 0}, {1, 3, 6, 2}, {1, 6, 4, 5}, {3, 4, 6, 7}, {1, 3, 4, 6} } };

                for (int tet=0; tet<5; ++tet)
                {
#ifdef VERBOSE
                    std::cout << "Tet #" << tet <<": ";
#endif
                    int parity = (i+j+k) % 2;
                    double tvpos[4][3];
                    double val[4];
                    for (int vtx=0; vtx<4; ++vtx)
                    {
                        tvpos[vtx][0] = i + pos[tets[parity][tet][vtx]][0];
                        tvpos[vtx][1] = j + pos[tets[parity][tet][vtx]][1];
                        tvpos[vtx][2] = k + pos[tets[parity][tet][vtx]][2];
                        val[vtx] = d[ptId[tets[parity][tet][vtx]]];
                    }
                    cs.addTetrahedron(tvpos, val);

                    PiecewisePolynomialFunction<double> csAreaPPF;
                    PiecewisePolynomialFunction<double> csVolPPF;
                    //ContourSpectrum<Identity>::calcSurfaceAreaAndVolumePFFForTet(tvpos, val, csAreaPPF, csVolPPF);
                    ContourSpectrum<RoundToNearestMultiple>::calcSurfaceAreaAndVolumePFFForTet(tvpos, val, csAreaPPF, csVolPPF, roundF);

                    unsigned int minValVtx = 0;
                    unsigned int maxValVtx = 0;
                    for (int vtx=0; vtx<4; ++vtx)
                    {
                        if (d[ptId[tets[parity][tet][vtx]]] < d[ptId[tets[parity][tet][minValVtx]]])
                            minValVtx = vtx;
                        if (d[ptId[tets[parity][tet][vtx]]] > d[ptId[tets[parity][tet][maxValVtx]]])
                            maxValVtx = vtx;

#ifdef VERBOSE
                        std::cout << int(d[ptId[tets[parity][tet][vtx]]]) << " ";
#endif
                    }
#ifdef VERBOSE
                    std::cout << "minValVtx=" << minValVtx << " maxValVtx=" << maxValVtx << " Path: ";
#endif
                    std::list<ctBranch*> path1, path2;
                    ctBranch* commonParent;
                    findPath(bm[ptId[tets[parity][tet][minValVtx]]],
                             bm[ptId[tets[parity][tet][maxValVtx]]],
                             path1, path2, commonParent);

                    //double currMinVal = d[ptId[tets[parity][tet][minValVtx]]];
                    //double currMaxVal = d[ptId[tets[parity][tet][maxValVtx]]];
                    double currMinVal;
                    if (path1.empty())
                    {
                        currMinVal = std::min(d[commonParent->extremum], d[commonParent->saddle]);
                    }
                    else
                    {
                        ctBranch *b= path1.front();
                        currMinVal = std::min(d[b->extremum], d[b->saddle]);
                    }
                    assert(currMinVal <= d[ptId[tets[parity][tet][minValVtx]]]);

                    double currMaxVal;
                    if (path2.empty())
                    {
                        currMaxVal = std::max(d[commonParent->extremum], d[commonParent->saddle]);
                    }
                    else
                    {
                        ctBranch *b= path2.back();
                        currMaxVal = std::max(d[b->extremum], d[b->saddle]);
                    }
                    assert(currMaxVal >= d[ptId[tets[parity][tet][maxValVtx]]]);

                    bool minInsideInterval = true;

                    for (std::list<ctBranch*>::iterator it = path1.begin(); it != path1.end(); ++it)
                    {
#ifdef VERBOSE
                        std::cout << "[ Branch = (" <<int(d[(*it)->saddle]) << ", " << int(d[(*it)->extremum]) << "), Range = (" <<  currMinVal << " , " << int(d[(*it)->saddle]) <<") ] ";
#endif
                        if (currMinVal < d[(*it)->saddle])
                        {
                            static_cast<BranchData*>((*it)->data)->areaPPF += csAreaPPF.restrictTo(currMinVal, d[(*it)->saddle]);

                            std::list<double> saddleValsOnPath;
                            ITERATE_OVER_CHILDREN(c, *it)
                            {
                                if (d[c->saddle] > currMinVal && d[c->saddle] > csVolPPF.domainMin())
                                {
                                    //assert(d[c->saddle] > csVolPPF.domainMin());
                                    assert(saddleValsOnPath.empty() || d[c->saddle] >= saddleValsOnPath.back());
                                    if (saddleValsOnPath.empty() || d[c->saddle] > saddleValsOnPath.back())
                                    {
                                        //std::cout << d[c->saddle] << " ";
                                        saddleValsOnPath.push_back(d[c->saddle]);
                                    }
                                }
                            }
                            //std::cout << std::endl;
                            assert(saddleValsOnPath.empty() || d[(*it)->saddle] >= saddleValsOnPath.back());
                            if (saddleValsOnPath.empty() || d[(*it)->saddle] > saddleValsOnPath.back())
                                saddleValsOnPath.push_back(d[(*it)->saddle]);
                            //saddleValsOnPath.sort();
                            //saddleValsOnPath.unique();
                            for (std::list<double>::iterator sit=saddleValsOnPath.begin();
                                    sit!=saddleValsOnPath.end(); ++sit)
                            {
                                assert(currMinVal != *sit);
                                PiecewisePolynomialFunction<double> addVolPPF = csVolPPF.restrictTo(currMinVal, *sit);
                                if (!addVolPPF.empty())
                                {
                                    if (!minInsideInterval) addVolPPF -= addVolPPF(addVolPPF.domainMin());
                                    static_cast<BranchData*>((*it)->data)->volumePPF += addVolPPF;
                                    minInsideInterval = false;
                                }
                                else
                                {
                                    std::cerr << "Warning: Empty vol ppf for restriction to interval [" << currMinVal << ", " << *sit << "]. volPPF is defined between " << csVolPPF.domainMin() << " and " << csVolPPF.domainMax() << std::endl;
                                }
                                currMinVal = *sit;
                            }
                        }
                    }
                    std::list<ctBranch*>::reverse_iterator rit;
                    for (rit = path2.rbegin(); rit != path2.rend(); ++rit)
                    {
#ifdef VERBOSE
                        std::cout << "[ Branch = (" <<int(d[(*rit)->saddle]) << ", " << int(d[(*rit)->extremum]) << "), Range = (" << int(d[(*rit)->saddle]) << " , " << currMaxVal << ") ] ";
#endif
                        if (currMaxVal > d[(*rit)->saddle])
                        {
                            static_cast<BranchData*>((*rit)->data)->areaPPF += csAreaPPF.restrictTo(d[(*rit)->saddle], currMaxVal);

			    std::list<double> saddleValsOnPath;
			    saddleValsOnPath.push_back(d[(*rit)->saddle]);

                            ITERATE_OVER_CHILDREN(c, *rit)
                            {
                                if (d[c->saddle] < currMaxVal && d[c->saddle] > csVolPPF.domainMin())
                                {
                                    //assert(d[c->saddle] > csVolPPF.domainMin());
                                    assert(saddleValsOnPath.empty() || d[c->saddle] >= saddleValsOnPath.back());
                                    if(saddleValsOnPath.empty() || d[c->saddle] > saddleValsOnPath.back())
                                    {
                                        //std::cout << d[c->saddle] << " ";
                                        saddleValsOnPath.push_back(d[c->saddle]);
                                    }
                                }
                            }
                            //std::cout << std::endl;
                            //saddleValsOnPath.sort();
                            //saddleValsOnPath.unique();
                            for (std::list<double>::reverse_iterator sit=saddleValsOnPath.rbegin();
                                    sit!=saddleValsOnPath.rend(); ++sit)
                            {
                                assert(*sit != currMaxVal);
                                PiecewisePolynomialFunction<double> addVolPPF = csVolPPF.restrictTo(*sit, currMaxVal);
                                if (!addVolPPF.empty())
                                {
                                    addVolPPF -= addVolPPF(addVolPPF.domainMin());
                                    static_cast<BranchData*>((*rit)->data)->volumePPF += addVolPPF;
                                }
                                else
                                {
                                    std::cerr << "Warnining: Empty vol ppf for restriction to interval [" << *sit << ", " << currMaxVal  << "]. volPPF is defined between " << csVolPPF.domainMin() << " and " << csVolPPF.domainMax() << std::endl;
                                }
                                currMaxVal = *sit;

                            }
                        }
                    }
#ifdef VERBOSE
                    std::cout << "[ Root Branch = ("<<int(d[commonParent->saddle]) << ", " << int(d[commonParent->extremum]) << "), Range = (" << currMinVal << ", " << currMaxVal << ") ]";
#endif
                    static_cast<BranchData*>(commonParent->data)->areaPPF += csAreaPPF.restrictTo(currMinVal, currMaxVal);

                    if (currMinVal != currMaxVal)
                    {
                        std::list<double> saddleValsOnPath;
                        //assert(currMinVal >= csVolPPF.domainMin());
                        ITERATE_OVER_CHILDREN(c, commonParent)
                        {
                            if (d[c->saddle] > currMinVal && d[c->saddle] > csVolPPF.domainMin() && d[c->saddle] < currMaxVal)
                            {
                                assert(saddleValsOnPath.empty() || d[c->saddle] >= saddleValsOnPath.back());
                                if(saddleValsOnPath.empty() || d[c->saddle] > saddleValsOnPath.back())
                                {
                                    //std::cout << d[c->saddle] << " ";
                                    saddleValsOnPath.push_back(d[c->saddle]);
                                }
                            }
                        }
                        //std::cout << std::endl;
                        assert(saddleValsOnPath.empty() || currMaxVal >= saddleValsOnPath.back());
                        if(saddleValsOnPath.empty() || currMaxVal > saddleValsOnPath.back())
                            saddleValsOnPath.push_back(currMaxVal);
                        //saddleValsOnPath.sort();
                        //saddleValsOnPath.unique();
                        for (std::list<double>::iterator sit=saddleValsOnPath.begin();
                                sit!=saddleValsOnPath.end(); ++sit)
                        {
                            assert(currMinVal != *sit);
                            PiecewisePolynomialFunction<double> addVolPPF = csVolPPF.restrictTo(currMinVal, *sit);
                            if (!addVolPPF.empty())
                            {
                                if (!minInsideInterval) addVolPPF -= addVolPPF(addVolPPF.domainMin());
                                static_cast<BranchData*>(commonParent->data)->volumePPF += addVolPPF;
                                minInsideInterval = false;
                            }
                            else
                            {
                                std::cerr << "Warnining: Empty vol ppf for restriction to interval [" << currMinVal << ", " << *sit  << "]. volPPF is defined between " << csVolPPF.domainMin() << " and " << csVolPPF.domainMax() << std::endl;
                            }
                            currMinVal = *sit;
                        }
                    }
#ifdef VERBOSE
                    std::cout << std::endl;
#endif
                }
            }
        std::cout << "." << std::flush;
    }
    std::cout << std::endl;
    std::cout << "Saving area & volume" << std::endl;

    cs.combine();
    cs.saveArea("total_area.dat", 200);
    cs.saveOutsideVolume("total_volume.dat", 200);
}

void saveForVisIt(std::ofstream &os, ctBranch *b, Data &d)
{
    typedef uint64_t IdxT;
    typedef double ValueT;

    IdxT extremum = b->extremum;
    os.write(reinterpret_cast<const char*>(&extremum), sizeof(extremum));
    ValueT extremumVal = d[extremum];
    os.write(reinterpret_cast<const char*>(&extremumVal), sizeof(extremumVal));
    IdxT saddle = b->saddle;
    os.write(reinterpret_cast<const char*>(&saddle), sizeof(extremum));
    ValueT saddleVal = d[saddle];
    os.write(reinterpret_cast<const char*>(&saddleVal), sizeof(saddleVal));
    IdxT volume = 0;
    os.write(reinterpret_cast<const char*>(&volume), sizeof(volume));
    IdxT id = static_cast<BranchData*>(b->data)->id;
    os.write(reinterpret_cast<const char*>(&id), sizeof(id));
    IdxT num = 0;
    ITERATE_OVER_CHILDREN(c, b)
        ++num;
    os.write(reinterpret_cast<const char*>(&num), sizeof(num));
    ITERATE_OVER_CHILDREN(c, b)
        saveForVisIt(os, c, d);
}

int main( int argc, char ** argv )
{
    int errflg = 0;
    int c;

    // Command line parameters
    char filename[1024] = "";
    char outfile[1024] = "";
    char branchmapfilename[1024] = "";
    char visittreefilename[1024] ="";
    bool computeContourSpectrum = false;

    char switches[256] = "i:o:b:v:c";

    while ((c = getopt(argc, argv, switches)) != EOF)
    {
        switch (c)
        {
            case 'i': {
                          strcpy(filename,optarg);
                          break;
                      }
            case 'o': {
                          strcpy(outfile,optarg);
                          break;
                      }
            case 'b': {
                          strcpy(branchmapfilename, optarg);
                          break;
                      }
            case 'v': {
                          strcpy(visittreefilename, optarg);
                          break;
                      }
            case 'c':
                      {
                          computeContourSpectrum = true;
                          break;
                      }
            case '?':
                      errflg++;
        }
    }

    if (errflg || filename[0] == '\0')
    {
        std::clog << "usage: " << argv[0] << " <flags> " << std::endl << std::endl;

        std::clog << "flags" << std::endl;
        std::clog << "\t -i < filename >  :  input filename" << std::endl;
        std::clog << "\t -o < filename >  :  output tree filename" << std::endl;
        std::clog << "\t -b < filename >  :  branch map filename" << std::endl;
        std::clog << "\t -v < filename >  :  VisIt compatible tree filename" << std::endl;
        std::clog << "\t -c : Compute contour spectrum" << std::endl;
        std::clog << std::endl;

        std::clog << "Filename must be of the form <name>.<i>x<j>x<k>.<type>" << std::endl;
        std::clog << "i,j,k are the dimensions. type is one of uint8, uint16, float or double." << std::endl << std::endl;
        std::clog << "eg, \" turtle.128x256x512.float \" is a file with 128 x 256 x 512 floats." << std::endl;
        std::clog << "i dimension varies the FASTEST (in C that looks like \"array[k][j][i]\")" << std::endl;
        std::clog << "Data is read directly into memory -- ENDIANESS IS YOUR RESPONSIBILITY." << std::endl;

        return(1);
    }

    if (computeContourSpectrum && outfile[0] == '\0')
    {
        std::clog << "Computing contour spectrum only makes sense when saving tree." << std::endl;
        return (2);
    }

    char prefix[1024];

    // Load data
    Data data;
    bool compress;
    if (!data.load(filename, prefix, &compress))
    {
        std::cerr << "Failed to load data" << std::endl;
        exit(1);
    }

    // Create mesh
    Mesh mesh(data);
    std::vector<size_t> totalOrder;
    mesh.createGraph(totalOrder); // This just sorts the vertices according to data.less()

    // Init libtourtre
    ctContext * ctx = ct_init(
            data.totalSize, //numVertices
            &(totalOrder.front()), //totalOrder. Take the address of the front of an stl vector, which is the same as a C array
            &value,
            &neighbors,
            &mesh //data for callbacks. The global functions less, value and neighbors are just wrappers which call mesh->getNeighbors, etc
            );

    // Create contour tree
    std::cout << "Computing contour tree and branch decomposition." << std::endl;
    ct_sweepAndMerge( ctx );
    ctBranch *root = ct_decompose(ctx);
    ctBranch **map = ct_branchMap(ctx);

    std::cout << "Simplifying tree and adding data structures." << std::endl;
    IdT lastUnsimplifiedId = addId(root);
    IdT *idBasedMap = createIdBasedBranchMap(map, data.totalSize);
    IdT *oldToNewIdMap = new IdT[lastUnsimplifiedId];
    IdT currId = 0;
    ctBranch * simplifiedRoot = ctBranch_copyRemovePerturbationAndAddBranchData(root, data, ctx, oldToNewIdMap, currId);
    applyIdTranslationToBranchMap(idBasedMap, data.totalSize, oldToNewIdMap);
    ctBranch** idToBranchPtrMap = createIdToBranchPtrMap(simplifiedRoot, currId);
    delete[] oldToNewIdMap;
    for (int idx=0; idx<data.totalSize; ++idx)
        map[idx] = idToBranchPtrMap[idBasedMap[idx]];

#ifndef NDEBUG
    checkTree(simplifiedRoot, data);
#endif

    outputTree(std::cout, simplifiedRoot, data);
    std::cout << "Calculating contour spectrum." << std::endl;
    if (computeContourSpectrum) calcContourSpectrum(data, map);
    delete[] idBasedMap;

    ct_cleanup( ctx );

    //output tree
    outputTree(std::cout, simplifiedRoot, data);
    std::cout << std::endl;

    std::ofstream os("tree_total_volume.dat");
    for (int i=0; i < 200; ++i)
    {
        double v = data.minValue + (i/200.)*(data.maxValue-data.minValue);
        os << v << " " << computeVolumeFromTree(simplifiedRoot, data, v) << std::endl;
    }
    os.close();
    os.open("tree_total_area.dat");
    for (int i=0; i < 200; ++i)
    {
        double v = data.minValue + (i/200.)*(data.maxValue-data.minValue);
        os << v << " " << computeAreaFromTree(simplifiedRoot, data, v) << std::endl;
    }

    if (branchmapfilename[0] != '\0')
    {
        FILE *bf = fopen(branchmapfilename, "w");
        for (size_t i = 0; i < data.totalSize; ++i)
            fwrite(static_cast<void*>(&(static_cast<BranchData*>(map[i]->data)->id)), sizeof(IdT), 1, bf);
        fclose(bf);
    }

    if (visittreefilename[0] != '\0')
    {
        std::ofstream os(visittreefilename);
        saveForVisIt(os, simplifiedRoot, data);
    }

    if (outfile[0] != '\0')
    {
        ctBranch_save(simplifiedRoot, outfile, computeContourSpectrum ? saveCS : 0);
    }
}
