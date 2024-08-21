#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define M 4 // the maximum number of indices a node can contain
#define m 2 // the minimum number of indices a node can contain
#define N 30 // the order of Hilbert curve

// Point
typedef struct point{ // point (x, y)
    long long int x; // x value
    long long int y; // y value
} Point;
// Point pointer
typedef Point* POINT;

// Rectangle 
typedef struct rectangle{
    POINT min; // bottom left vertex
    POINT max; // upper right vertex
} Rectangle;
// Rectangle pointer
typedef Rectangle* RECTANGLE;

// Index Pointer
typedef struct index* INDEX;
// Node Pointer
typedef struct node* NODE;

// Node
typedef struct node{
    NODE parent; // parent node of this node
    INDEX I[M]; // the array of indices stored in this node
    int size; // how many indices are actually stored inside the array
} Node;

// Index 
typedef struct index{
    RECTANGLE rect; // the MBR of all rectangles present in the child node of this index
    long long int lhv; // the Largest Hilbert Value among child indices
    NODE child; // a pointer to the child node of this index
} Index;

// Tree
typedef struct tree{
    NODE root;
} Tree;
// Tree Pointer
typedef struct tree* TREE;

// constructor method for a point
POINT newPoint(long long int x, long long int y){
    POINT newPoint = (POINT) malloc(sizeof(Point));
    newPoint->x = x;
    newPoint->y = y;
    return newPoint;
}

// constructor method for a index
INDEX newIndex(RECTANGLE rect, long long int lhv, NODE child){
    INDEX newIndex = (INDEX) malloc(sizeof(Index));
    newIndex->rect = rect;
    newIndex->lhv = lhv;
    newIndex->child = child;
    return newIndex;
}

// constructor method for a rectangle
RECTANGLE newRectangle(POINT min, POINT max){
    RECTANGLE newRect = (RECTANGLE) malloc(sizeof(Rectangle));
    newRect->min = min;
    newRect->max = max;
    return newRect;
}

// constructor method for a new node
NODE newNode(){
    NODE newNode = (NODE) malloc(sizeof(Node));
    for (int i = 0; i < M; i++) newNode->I[i] = NULL;
    newNode->size = 0;
    newNode->parent = NULL;
    return newNode;
}

// constructor method for a new tree
TREE newTree(){
    NODE new = newNode();
    TREE newTree = (TREE) malloc(sizeof(Tree));
    newTree->root = new;
    return newTree;
}

char newQuadType(char oldType, long long int quad){
    if(oldType=='A')
    {
        if(quad==0) return 'D';
        if(quad==1) return 'A';
        if(quad==2) return 'A';
        if(quad==3) return 'B';
    }
    if(oldType=='B')
    {
        if(quad==0) return 'B';
        if(quad==1) return 'B';
        if(quad==2) return 'C';
        if(quad==3) return 'A';
    }
    if(oldType=='C')
    {
        if(quad==0) return 'C';
        if(quad==1) return 'D';
        if(quad==2) return 'B';
        if(quad==3) return 'C';
    }
    if(oldType=='D')
    {
        if(quad==0) return 'A';
        if(quad==1) return 'C';
        if(quad==2) return 'D';
        if(quad==3) return 'D';
    }
    return 'A';
}

long long int quadValue(char oldType, long long int quad){
    if(oldType=='A')
    {
        if(quad==0) return 0;
        if(quad==1) return 1;
        if(quad==2) return 2;
        if(quad==3) return 3;
    }
    if(oldType=='B')
    {
        if(quad==0) return 2;
        if(quad==1) return 1;
        if(quad==2) return 0;
        if(quad==3) return 3;
    }
    if(oldType=='C')
    {
        if(quad==0) return 2;
        if(quad==1) return 3;
        if(quad==2) return 0;
        if(quad==3) return 1;
    }
    if(oldType=='D')
    {
        if(quad==0) return 0;
        if(quad==1) return 3;
        if(quad==2) return 2;
        if(quad==3) return 1;
    }
    return 0;
}

void hilbert(long long int n, long long int x, long long int y, long long int* res, char quadType){
    if(n==1) return;
    long long int p=n/2;
    long long int quad;
    if(x<p && y<p)
    {
        quad=0;
    }
    else if(x<p && y>=p)
    {
        quad=1;
        y-=p;
    }
    else if(x>=p && y>=p)
    {
        quad=2;
        x-=p;
        y-=p;
    }
    else
    {
        quad=3;
        x-=p;
    }
    long long int quadVal=quadValue(quadType,quad);
    *res=((*res)*pow(2, 2))+(quadVal);
    char newQuad=newQuadType(quadType, quad);
    hilbert(n/2,x,y,res, newQuad);
}

long long int hilbertValue(POINT p){
    long long int res=0;
    hilbert((long long int)pow(2.0, N), p->x, p->y, &res, 'A');
    return res;
}

// checks if a node is a leaf node
int isLeaf(NODE n){
    return n->I[0]->child == NULL;
}

// checks if thwo rectangles intersect in 2D space
int checkIntersection(RECTANGLE r1, RECTANGLE r2){
    return (r2->min->x <= r1->max->x) && (r2->max->x >= r1->min->x) && (r2->min->y <= r1->max->y) && (r2->max->y >= r1->min->y);
}

// returns the leaf node in which to place a new rectangle
NODE ChooseLeaf(NODE n, RECTANGLE r, long long int h){
    //n is the root node at which to start
    if(isLeaf(n)){
        return n;
    }
    for(int i=0; i<n->size; i++){
        if(n->I[i]->lhv > h){
            return ChooseLeaf(n->I[i]->child, r, h);
        }
    }
    return ChooseLeaf(n->I[n->size-1]->child, r, h); // 
}

// returns the cooperating sibling of a node
NODE GetCooperatingSibling(NODE n){
    NODE parent = n->parent;
    NODE sibling = NULL;
    // n can only have siblings if its parent exists and has a size of at least 2
    if(parent && (parent->size > 1)){ 
        int i=0; 
        while(parent->I[i]->child != n) i++; // get the index of the MBR of n in its parent

        if(i==0) sibling = parent->I[i+1]->child; // if it is the first in line, we take the node to its right
        else if(i==parent->size-1) sibling = parent->I[i-1]->child; // if it is the last in line, we take the immediate left
        else{ // otherwise we check for both left and right siblings
            if(parent->I[i-1]->child->size < M) sibling = parent->I[i-1]->child; // incase the left sibling is not full
            else sibling = parent->I[i+1]->child; // in case left is full, we take the right one regardless it is full or not
        }
    }
    return sibling;
}

// handles insertion in an overflowing node by
// either splitting it into two or taking the 
// help of a cooperating sibling node to split
NODE HandleOverflow(TREE t, NODE n, INDEX ind){
    NODE root = t->root;

    // get the sibling node 
    NODE sibling = GetCooperatingSibling(n);
    // now we calculate the size of the set of indices 
    // present in n, its sibling and the new index
    int setSize = n->size + 1; 
    if(sibling) setSize += sibling->size;

    INDEX set[setSize];
    // add indices of n and sibling in sorted order
    int j = 0;

    if(sibling){
        if(n->I[0]->lhv < sibling->I[0]->lhv){    // if indices in n are smaller, append them first
            for(int i = 0; i<n->size; i++)
                set[j++] = n->I[i];
            for(int i = 0; i<sibling->size; i++)
                set[j++] = sibling->I[i];
        }
        else{                                     // otherwise, append the indices of the sibling first
            for(int i = 0; i<sibling->size; i++)
                set[j++] = sibling->I[i];
            for(int i = 0; i<n->size; i++)
                set[j++] = n->I[i];
        }
    }
    else{
        for(int i = 0; i<n->size; i++)           // if there is no sibling, just append the indices of n
                set[j++] = n->I[i];
    }
    // insert the new rectangle in sorted order
    if(ind->lhv > set[setSize - 2]->lhv)
        j = setSize-1;
    else{
        for(int i = 0; i<setSize-1; i++){
            j = i;
            if (set[i]->lhv > ind->lhv){  
                break;  // got the index
            }
        }
    }
    for (int i = setSize-2; i>=j; i--) //last index of set array is setSize-1
        set[i+1] = set[i];    // shift the elements(rectangles)
    set[j] = ind;             // insert the new rectangle
    
    // now we have our sorted set of rectangles.
    // Now, if there is space in the sibling, distribute all rectangles evenly in both
    if(sibling && (sibling->size < M)){ // now total elements can be adjusted in the two siblings
        int k = setSize/2;
        
        if(n->I[0]->lhv < sibling->I[0]->lhv){ // sibling is greater in terms of Hilbert Value
            int i1 = 0, i2 = 0, i3 = 0;
            while(i1 < k) 
                n->I[i1++] = set[i2++];  // set the node first
            while(i1 < n->size)
                n->I[i1++] = NULL;
            while(i2 < setSize)    // then set the sibling
                sibling->I[i3++] = set[i2++];
            while(i3 < sibling->size)
                sibling->I[i3++] = NULL;
            n->size = k; 
            sibling->size = setSize - k;
        }
        else{
            int i1 = 0, i2 = 0, i3 = 0;
            while(i3 < k)
                sibling->I[i3++] = set[i2++];  // set the sibling first
            while(i3 < sibling->size)
                sibling->I[i3++] = NULL;
            while(i2 < setSize)
                n->I[i1++] = set[i2++];     // then set the node 
            while(i1 < n->size)
                n->I[i1++] = NULL;
            sibling->size = k;
            n->size = setSize - k;
        }

        return NULL;
    }
    else if(sibling){ // if sibling is present and it is full then we will distribute all indeices in 3 nodes
        NODE NN = newNode(); // new node
        int k = setSize/3;
        if(n->I[0]->lhv < sibling->I[0]->lhv){ // so sibling was on right 
            for(int i = 0; i < n->size; i++) // removing all indices from n
                n->I[i] = NULL;
            for(int i = 0; i < sibling->size; i++) // removing all indices from siblings
                sibling->I[i] = NULL;
            int j = 0;
            for(int i = 0; i < k ; i++) // putting first k indices in n (lowest) 
                n->I[i] = set[j++];
            for(int i = 0; i < k ; i++) // putting next k indices in siblings
                sibling->I[i] = set[j++];
            for(int i = 0; i < setSize - 2*k; i++) // putting next k in new n
                NN->I[i] = set[j++];
        }
        else{ // sibling on left
            for(int i = 0; i < n->size; i++) // removing all indices from n
                n->I[i] = NULL;
            for(int i = 0; i < sibling->size; i++) // removing all indices from siblings
                sibling->I[i] = NULL;
            int j = 0;
            for(int i = 0; i < k ; i++) // putting next k indices in siblings (lowest)
                sibling->I[i] = set[j++];
            for(int i = 0; i < k ; i++)// putting first k indices in n
                n->I[i] = set[j++];
            for(int i = 0; i < setSize - 2*k ; i++) // putting next k in new n
                NN->I[i] = set[j++];
        }
        // updating all the sizes
        n->size = k;
        sibling->size=k;
        NN->size=setSize - 2*k;
        
        return NN;    
    }
    else{ // when there is no sibling. This case should only arise when n is the root node
        NODE NN = newNode(); // create a new node to distribute indices evenly
        int k = setSize/2;   // k is the number of indices we will put in n, the original node
        
        for(int i = 0; i < n->size; i++) 
            n->I[i] = NULL;              // make all indices in n NULL
        n->size = 0;
        
        int j = 0;
        for(int i = 0; i < k; i++){
            n->I[i] = set[j++];     // add the first k indices of the set to n
            n->size++;
        }
        for(int i = 0; i < setSize - k; i++){
            NN->I[i] = set[j++];           // add the rest of them to new node NN
            NN->size++;
        }
        
        return NN;
    }
}

// returns the MBR of a node
INDEX MBR(NODE n){
    INDEX ind = NULL;
    if(n){
        POINT min = newPoint(0, 0), max = newPoint(0, 0);
        long long int lhv = n->I[0]->lhv;
        
        min->x = n->I[0]->rect->min->x;
        min->y = n->I[0]->rect->min->y;
        max->x = n->I[0]->rect->max->x;
        max->y = n->I[0]->rect->max->y;

        // calculate the min and max points and lhv of the MBR of n
        for(int i = 1; i < n->size; i++){
            if(n->I[i]->rect->min->x < min->x)
                min->x = n->I[i]->rect->min->x; // x value for min will be the minimum x for all rectangles in n
            if(n->I[i]->rect->min->y < min->y)
                min->y = n->I[i]->rect->min->y; // y value for min will be the minimum y for all rectangles in n
            if(n->I[i]->rect->max->x > max->x)
                max->x = n->I[i]->rect->max->x; // x value for max will be the maximum x for all rectangles in n
            if(n->I[i]->rect->max->y > max->y)
                max->y = n->I[i]->rect->max->y; // y value for max will be the maximum y for all rectangles in n
            if(n->I[i]->lhv > lhv)
                lhv = n->I[i]->lhv;             // lhv will be the maximum lhv for all rectangles in n
        }
        RECTANGLE mbr = newRectangle(min, max); // creating the MBR
        ind = newIndex(mbr, lhv, n); // creating the index that stores the MBR
    }
    return ind;
}

// goes up the tree and fixes the
// MBR and LHV values for nodes.
// Returns a new node if root is split 
NODE AdjustTree(TREE t, NODE arr[]){
    NODE root = t->root;
    /*
        arr[0] -> node being updated
        arr[1] -> sibling node, if any (can be NULL)
        arr[2] -> new node created during split, if any (can be NULL)
    */
    if(!arr[0]->parent) return arr[2]; // if reached root level, stop and return the new node if root was split
    NODE Np = arr[0]->parent;
    
    for(int i = 0; i < Np->size; i++){
        if(Np->I[i]->child == arr[0] || Np->I[i]->child == arr[1]) 
            Np->I[i] = MBR(Np->I[i]->child); // fix the MBRs present in the parent node
    }

    //Sorting the parent node in an adaptive way
    for(int i = 1; i < Np->size; i++){
        INDEX val = Np->I[i];
        int j;
        for(j = i - 1;j >= 0 && Np->I[j]->lhv > val->lhv;j--){
            Np->I[j+1] = Np->I[j];
        }
        Np->I[j+1] = val;
    }

    NODE parents[3]; 
    // parents is  an array that contains the parent node 
    // of the node being updated, its cooperating sibling,
    // if any, and the new node that might be formed when 
    // we try to insert into the parent node
    parents[0] = Np;
    parents[1] = GetCooperatingSibling(Np);
    parents[2] = NULL;
    
    // if there was a split when trying to insert into 
    // the node, we try to insert the MBR of the newly
    // created node into the parent node Np
    if(arr[2]){
        INDEX NN = MBR(arr[2]);
        
        if(Np->size < M){ // if there is space in Np, we insert
            for(int i = Np->size-1 ; i >= 0; i--){
                if(Np->I[i]->lhv > NN->lhv)
                    Np->I[i+1] = Np->I[i];
                else{
                    Np->I[i+1] = NN;
                    break;
                } 
            }
            Np->size++;
            arr[2]->parent = Np; // setting the parent of the new node
        }
        else{ // otherwise, call HandleOverflow
            parents[2] = HandleOverflow(t, Np, NN);

            // check the parent, its sibling, and the new parent node 
            // formed in HandleOverflow for the nodes present in arr.
            // Set the parent node of those nodes accordingly. 
            
            for(int k = 0; k < 3; k++){ // for all nodes in arr
                for(int i = 0; i < 3; i++){ // check all nodes in parents
                    if(parents[i]){
                        for (int j = 0; j < parents[i]->size; j++) // and all indices in those nodes
                            parents[i]->I[j]->child->parent = parents[i]; // assign the parent for all nodes
                    }
                }
            }
        }
    }

    return AdjustTree(t, parents);
}

// creates a new root node and increases
// the height of the tree if the root was
// split in two
void JoinRoot(TREE t, NODE left){
    NODE root = t->root; 
    NODE new = newNode();
    if(root->I[root->size - 1]->lhv < left->I[left->size - 1]->lhv){
        new->I[new->size++] = MBR(root);
        new->I[new->size++] = MBR(left);
    }
    else{
        new->I[new->size++] = MBR(left);
        new->I[new->size++] = MBR(root);
    }
    t->root = new;
    root->parent = new;
    left->parent = new;
}

// inserts a rectangle into a Hilbert R-Tree
void Insert(TREE t, RECTANGLE r){   
    NODE root = t->root;
    // calculating the center point of the rectangle
    long long int xmid=(r->min->x+r->max->x)/2;
    long long int ymid=(r->min->y+r->max->y)/2;
    POINT pmid=newPoint(xmid, ymid); // creating a point with those coordinates

    INDEX ind = newIndex(r, hilbertValue(pmid), NULL); // creating an index record having the given rectangle
    
    if(t->root->I[0] == NULL){
        t->root->I[0]=ind;
        t->root->size++;
        return;
    }
    // this selects the leaf node in which to place our rectangle
    NODE leaf = ChooseLeaf(root, r, hilbertValue(pmid)); 
    
    NODE arr[3];
    arr[0] = leaf;
    arr[1] = GetCooperatingSibling(leaf);
    arr[2] = NULL;

    if(leaf->size<M){ // if there is space in this node, insert the rectangle
        for(int i = leaf->size++; i>=0; i--){
            if(i > 0 && leaf->I[i-1]->lhv > ind->lhv) leaf->I[i] = leaf->I[i-1]; // inserting the rectangle in the
            else{                                                                // node in a sorted order based on
                leaf->I[i]=ind;                                                  // the hilbert values of all the 
                break;                                                           // rectangles, akin to the
            }                                                                    // insertion sort algorithm 
        }
    }
    else{ // else, call the HandleOverflow method
        arr[2] = HandleOverflow(t, leaf, ind);
    }
    
    NODE Left = AdjustTree(t, arr);
    if(Left)
        JoinRoot(t, Left);
}

// recursively searches a Hilbert R-Tree 
// for rectangles that intersect with the
// query rectangle, and prints out the data
// points that are bounded by it
void SearchInTree(NODE root, RECTANGLE r){
    for(int i=0; i<root->size; i++){
        if(checkIntersection(root->I[i]->rect, r)){
            if(isLeaf(root))
                printf("(%lld, %lld)-->", root->I[i]->rect->min->x, root->I[i]->rect->min->y);
            else
                SearchInTree(root->I[i]->child, r);
        }
    }
    printf("\n");  // print new line whenever one node is traversed
    return;
}

// wrapper for SearchInTree
void Search(TREE t, RECTANGLE r){
    SearchInTree(t->root, r);
}

// traverses the tree in a pre-order fashion
// and prints the node visited
void Traversal(NODE n){
    if(n){
        if(isLeaf(n)){
            printf("Leaf Node\n");
            for(int i = 0; i < n->size; i++){
                INDEX ind = n->I[i];
                printf("2D Data Point: (%lld, %lld)| %lld\n", ind->rect->min->x, ind->rect->min->y, ind->lhv);
            }
        }
        else{
            printf("Internal Node\n");
            for(int i = 0; i < n->size; i++){
                INDEX ind = n->I[i];
                printf("Bottom-Left Point: (%lld, %lld) ", ind->rect->min->x, ind->rect->min->y);
                printf("Top-Right Point: (%lld, %lld)| %lld\n", ind->rect->max->x, ind->rect->max->y, ind->lhv);
            }
        }
        printf("\n");
        for(int i = 0; i < n->size; i++) Traversal(n->I[i]->child);
    }
}

// wrapper for traversal
void PreOrderTraversal(TREE t){
    Traversal(t->root);
}

// driver code
int main(){

    char* filename = "DataSet.txt";               // change filename to test on a different file
    
    TREE t = newTree();                           // a new tree

    FILE* fp = fopen(filename, "r");              
    long long int x,y;
    while(fscanf(fp, "%lld %lld", &x, &y) == 2){  // read file for data while it is present
        POINT min=newPoint(x,y);                  // construct a point
        RECTANGLE r=newRectangle(min,min);        // a new rectangle (it is practically a point)
        Insert(t,r);                              // insert rectangle in the tree
    }
    fclose(fp);

    // sample call for traversal
    //
    PreOrderTraversal(t);
    // 

    // sample call for searching
    //
    // POINT qmin = newPoint(0,6), qmax = newPoint(10,18);
    // RECTANGLE q = newRectangle(qmin, qmax);
    // Search(t, q);
    //

    return 0;
}