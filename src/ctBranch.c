/*
Copyright (c) 2006, Scott E. Dillard
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
 are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation
    and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


#include "tourtre.h"
#include "ctBranch.h"
#include "ctContext.h"

ctBranch * ctBranch_new( size_t e, size_t s, ctContext * ctx )
{
    ctBranch * b = (*(ctx->branchAlloc))(ctx->cbData);
    b->extremum = e;
    b->saddle = s;
    b->parent = NULL;
    b->children = ctBranchList_init();
    b->nextChild = b->prevChild = NULL;
    return b;
}

/* recursive delete */
void ctBranch_delete( ctBranch * self, ctContext * ctx )
{ 
    ctBranch * c;
    for ( c = self->children.head; c != NULL; c = c->nextChild ) {
        ctBranch_delete( c, ctx );
    }
    (*(ctx->branchFree))(self,ctx->cbData);
}

ctBranchList ctBranchList_init()
{
    ctBranchList bl = { NULL };
    return bl;
}


static int 
compareSaddles(size_t a, size_t b, ctContext * ctx)
{
    return (*(ctx->value))(a,ctx->cbData) < (*(ctx->value))(b,ctx->cbData);
}

void ctBranchList_add(ctBranchList * self, ctBranch * c, ctContext * ctx)
{
    if (!self->head) { /* front of list */
        self->head = c;
        c->prevChild = NULL;
        c->nextChild = NULL;
    } else {
        ctBranch * i = self->head;
        while( compareSaddles(i->saddle,c->saddle,ctx) && i->nextChild ) 
            i = i->nextChild;
        
        if ( compareSaddles(i->saddle,c->saddle,ctx) ) { /* end of list */
            c->nextChild = NULL;
            c->prevChild = i;
            i->nextChild = c;
        } else { /* middle of list */
            c->nextChild = i;
            c->prevChild = i->prevChild;
            if (i->prevChild) i->prevChild->nextChild = c;
            else self->head = c;
            i->prevChild = c;
        }
    }
}
		

void ctBranchList_remove(ctBranchList * self, ctBranch * c)
{
    if (self->head == c) self->head = c->nextChild;
    if (c->nextChild) c->nextChild->prevChild = c->prevChild;
    if (c->prevChild) c->prevChild->nextChild = c->nextChild;
    c->nextChild = c->prevChild = NULL;
}


void ctBranchList_merge( ctBranchList * self, ctBranchList * other, ctContext * ctx )
{
    if (!other->head) return;
    if (!self->head) {
        self->head = other->head;
        return;
    }
    
    /* merge sort */
    {
        ctBranch *i = self->head;
        ctBranch *o = other->head;
        ctBranch *li = 0, *lo = 0;
        while(i && o) {
            if ( compareSaddles(o->saddle,i->saddle,ctx) ) {
                ctBranch * next;
                lo = o;
                /* move marker along */
                next = o->nextChild;
                /* insert oo into list before i */
                o->nextChild = i;
                o->prevChild = i->prevChild;
                if (i->prevChild) i->prevChild->nextChild = o;
                else self->head = o;
                i->prevChild = o;
                o = next;
            } else {
                li = i;
                i = i->nextChild;
            }
        }
        
        if (!i) {
            li->nextChild = o;
            o->prevChild = li;
        } else if (!o) {
            lo->nextChild = i;
            i->prevChild = lo;
        }
    }
}

void ctBranch_saveFile(ctBranch *b, FILE* f,
        void (*saveUserDataFct)(void *u, FILE* f))
{
    ctBranch *c;
    fwrite("(", sizeof(char), 1, f);
    fwrite(&(b->extremum), sizeof(size_t), 1, f);
    fwrite(&(b->saddle), sizeof(size_t), 1, f);
    if (saveUserDataFct) (*saveUserDataFct)(b->data, f);
    for (c = b->children.head; c != NULL; c = c->nextChild)
    {
        ctBranch_saveFile(c, f, saveUserDataFct);
    }
    fwrite(")", sizeof(char), 1, f);
}

void ctBranch_save(ctBranch *b, const char *filename,
        void (*saveUserDataFct)(void *u, FILE* f))
{
    FILE *f = fopen(filename, "w");
    if (f)
        ctBranch_saveFile(b, f, saveUserDataFct);
    else
        fprintf(stderr, "Error: Could not open file %s for writing.", filename);
}

ctBranch* ctBranch_loadFile(ctContext *ctx, FILE *f,
        void *(*loadUserDataFct)(FILE* f))
{
    char buff;
    size_t extremum, saddle;
    ctBranch *r;
    ctBranch *c;

    {
        fread(&extremum, sizeof(size_t), 1, f);
        fread(&saddle, sizeof(size_t), 1, f);
        r = ctBranch_new(extremum, saddle, ctx);
        if (loadUserDataFct)
            r->data = (*loadUserDataFct)(f);
        else
            r->data = 0;
        fread(&buff, sizeof(char), 1, f);
        while (buff == '(')
        {
            c = ctBranch_loadFile(ctx, f, loadUserDataFct);
            if (!c)
            {
                ctBranch_delete(r, ctx);
                return 0;
            }
            else
            {
                ctBranchList_add(&(r->children), c, ctx);
		c->parent = r;
            }
            fread(&buff, sizeof(char), 1, f);
        }
        if (buff != ')')
        {
            fprintf(stderr, "Error: Corrupt file. Expected branch start or end marker.\n");
            ctBranch_delete(r, ctx);
            return 0;
        }

        return r;
    }
}

ctBranch* ctBranch_load(ctContext *ctx, const char* filename,
        void *(*loadUserDataFct)(FILE* f))
{
    char buff;
    FILE *f;
    f = fopen(filename, "r");
    if (f)
    {
        fread(&buff, sizeof(char), 1, f);
        if (buff == '(')
        {
            return ctBranch_loadFile(ctx, f, loadUserDataFct);
        }
        else
        {
            fprintf(stderr, "Error: Corrupt file. Expected branch start marker\n");
            return 0;
        }
    }
    else
    {
        fprintf(stderr, "Error: Could not open file %s for reading.", filename);
        return 0;
    }
}
