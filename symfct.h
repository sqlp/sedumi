#if !defined(SYMFCT)
#define SYMFCT

typedef mwSignedIndex integer;  /* removed "long" */

int sfinit_(integer *neqns, integer *nnza, integer *xadj, integer *adjncy,
            integer *perm, integer *invp, integer *colcnt, integer *nnzl,
            integer *nsub, integer *nsuper, integer *snode, integer *xsuper,
            integer *iwsiz, integer *iwork, integer *iflag);

int symfct_(integer *neqns, integer *adjlen, integer *xadj, integer *adjncy,
            integer *perm, integer *invp, integer *colcnt, integer *nsuper,
            integer *xsuper, integer *snode, integer *nofsub, integer *xlindx,
            integer *lindx, integer *xlnz, integer *iwsiz, integer *iwork,
            integer *flag__);

#endif
