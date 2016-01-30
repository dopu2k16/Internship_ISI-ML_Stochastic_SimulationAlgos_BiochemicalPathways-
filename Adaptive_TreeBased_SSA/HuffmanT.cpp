#include <iostream>
#include <queue>
#include <string.h>
#include <vector>
#include<list>
#include <algorithm>
#include<new>
#define m 5
//#define size 2*m-1
using namespace std;

int Table[5]; //Rxn propensity vector

struct Node
{
    int prop;  //propensity value
    int rxno;  //reaction number given propensity value
    Node    * left;  //left & right pointer
    Node    * right;
    Node():left(NULL), right(NULL) , prop(0),rxno(0) {}
};


class Cmp   //declaring class for comparing the pointers to create a min heap
{
public:
    bool operator() (const Node * a, const Node * b) const
    {
        return  a->prop > b->prop;
    }
};

priority_queue<Node* , vector<Node*> , Cmp> myQueue; //for min heap & Queue operation

void BuildTree() //Building Huffman Tree
{

    for (int i = 0; i < m; ++ i)

        {
            Node * node = new Node(); //creating a new node of rxn propensity value
            node->prop = Table[i];
            node->left=nullptr;
            node->right=nullptr;
            node->rxno=i+1;


            myQueue.push(node);  //pushing the contents into priority queue which is min -heap
        }
    cout << "qUEQUE Top  "<<myQueue.top()->prop<<endl;


    while (myQueue.size() > 1)  //Building Huffman Tree
    {
        Node * f = myQueue.top();
        myQueue.pop();
        Node * s = myQueue.top();
        myQueue.pop();
        Node * tmp = new Node();
        tmp->prop = f->prop+ s->prop;
        cout<<tmp->prop<<endl;
        tmp->left = f;
        tmp->right = s;

       cout<<"Hello 2"<<endl;

        myQueue.push(tmp);
    }

cout << "qUEQUE Top"<<myQueue.top()->prop<<endl;
}

unsigned searchT( Node* huff,int val)  //search function for Huffman tree
{

//size=2*m-1;
   if(huff->left==nullptr && huff->right==nullptr)
        return huff->rxno ;

        if(huff->left->prop >= val)
          return   searchT(huff->left,val);
        else
        {

            //pos=2*start+2;
            val=huff->prop-val;
            return searchT(huff->right,val);

        }

    }



int main()
{
    //int N;
    memset(Table , 0 , sizeof(Table));
    //cin >> N;
    for (int i = 0; i < m; ++i)
    {
        cin>>Table[i];
    }

    BuildTree();
    //searchT()

 int val=185;
  unsigned sel=searchT(myQueue.top(),val);
  cout<<"SEL "<<sel<<endl;


    return 0;
}

