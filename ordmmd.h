#if !defined(ORDMMD)
#define ORDMMD

typedef mwSignedIndex integer;  /* removed "long" */

int ordmmd_(integer *neqns, integer *xadj, integer *adjncy, integer *invp,
            integer *perm, integer *iwsiz, integer *iwork, integer *nofsub,
            integer *iflag);

#endif
