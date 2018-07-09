#include<iostream>

# define Nkey 8

/* Implements a basic splay bst -- currently deletion doesn't splay  */
/* Splaying promotes the most recently added/queried node to root    */

typedef struct Bnode
{
int key[Nkey];         
int dof;
int ndof;
struct Bnode *left;
struct Bnode *right;
} bnode;

enum BstBool {BstLT,BstGT,BstEQ};

int bstadd(bnode *(*root), int key[], int ndof); /* return 0/+ndof */
int bstdat(bnode *(*root), int key[]);   /* Returns dof or -1   */

int  bstset(bnode *root);                /* assigns the node numbering */
void bstset0(bnode *root, int &dof);     /* helper for bstset()        */

int* bstsort(int key[]);                 /* sort descending MODIFIES key[] */
void bstdel(bnode *(*root), int key[]);  /* Deletes an entry     */
void bstfree(bnode *root);

bnode *bstcutmin(bnode *(*bst));         /* cuts and returns min */
int    bstcmp(int key0[], int key1[]);   /* compares two bst keys*/

int bstcmp(int k0[], int k1[])
{
int i;

for(i = 0; i < Nkey; i++)
{
if(k0[i] < k1[i]) return(BstLT);
if(k0[i] > k1[i]) return(BstGT);
}

return BstEQ;
}

int* bstsort(int pp[])   /* sort descending */
{
int i,j,gap,temp, n=Nkey;

for(gap = n/2; gap > 0; gap /=2)
for(i = gap; i < n; i++)
for(j = i-gap; j >= 0 && pp[j] < pp[j+gap]; j -= gap)
{
temp = pp[j];
pp[j] = pp[j+gap];
pp[j+gap] = temp;
}

return pp;
}

int bstadd(bnode *(*bst), int key[], int ndof)
{
int i,j;
bnode *root = *bst, *v;

if(root != NULL)
   {
   switch (bstcmp(key,root->key))
     {
     case BstLT:
       {
       i = bstadd(&(root->left), key, ndof);
	  
       v = root->left;		    /* promote left node	*/
       root->left = v->right;
       v->right = root;
	  
       *bst = v;
	  
       return i;
       }
	  
     case BstGT:
       {
       i = bstadd(&(root->right), key, ndof);
	
       v = root->right;		    /* promote right node	*/
       root->right = v->left;
       v->left = root;
	
       *bst = v;
	
       return i;
       }
	
     case BstEQ:
       {	  
       if(ndof == root->ndof) return(0);

       cout << "Error with bstadd(): "
	    << "nodf = " << root->ndof << ", requested ndof = "
            << ndof << endl;

       throw("bst error");
       }
     }
   }
 else
   {
   root = (bnode *) new bnode;
    
   for(j=0; j<Nkey; j++) root->key[j] = key[j];

   root->ndof = ndof;
   root->left = NULL;
   root->right= NULL;
    
   *bst = root;
    
   return ndof;
   }

 throw("bstadd(): error");
}

void bstdel(bnode *(*bst), int key[])  /* doesn't splay */
{
 bnode *root = *bst, *v;
 
 switch (bstcmp(key, root->key))
   {		 
   case BstLT: 
     {
     bstdel(&(root->left ), key);
	
     return;
     }

   case BstGT: 
     {
     bstdel(&(root->right), key);

     return;
     }

   case BstEQ:
     {
     if((root->right == NULL) 
	&&  (root->left  == NULL))*bst = NULL;
     else if(root->right == NULL) *bst = root->left;
     else if(root->left  == NULL) *bst = root->right;
     else
       {
       v = bstcutmin(&(root->right));
       *bst = v;
   
       v->left  = root->left;
       v->right = root->right;
       }

     delete root;
     return;
     }
   }
}

bnode *bstcutmin(bnode *(*bst))
{
 bnode *root = *bst;

 if(root->left == NULL)
   {
   *bst = root->right;
   
   return root;
   }
 else
   {
   return bstcutmin(&(root->left));
   }
}

int bstdat(bnode *(*bst), int key[])  /* returns the data for the key */
{
 int i;
 bnode *root = *bst, *v;

 if(root != NULL)
   {
   switch (bstcmp(key, root->key))
     { 
     case BstLT:
       {
       i = bstdat(&(root->left), key);

       if(i == -1) return i;  /* key doesn't exist   */

       v = root->left;		 /* promote left node	*/
       root->left = v->right;
       v->right = root;
	
       *bst = v;
	
       return i;
       }

     case BstGT:
       {
       i = bstdat(&(root->right), key);
	
       if(i == -1) return i;

       v = root->right;		 /* promote right node	*/
       root->right = v->left;
       v->left = root;
	
       *bst = v;

       return i;
       }

     case BstEQ:
       {
       return root->dof;
       }
     }
   }
 else
   {
   return -1;    /* key doesn't exist */
   }

 throw("bstdat(): error");
}

void bstscn(bnode *root, void f(void *, int key[], int dof), void* farg)
{
 if(root == NULL) return;

 bstscn(root->left, f, farg);

 f(farg, root->key, root->dof);

 bstscn(root->right, f, farg);

 return;
}

void bstset0(bnode *root, int &dof)
{
  if(root == NULL) return;

  bstset0(root->left, dof);

  root->dof = dof;
  dof += root->ndof;

  bstset0(root->right, dof);

  return;
}

int bstset(bnode *root)
{
  int dof = 0;

  bstset0(root, dof);

  return(dof);
}

void bstfree(bnode *root)
{
 if(root != NULL)
   {
   bstfree(root->left);
   bstfree(root->right);

   delete root;
   }

 return;
}
