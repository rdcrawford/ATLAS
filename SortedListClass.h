#include <iostream>
#include "LinkedNodeClass.h"
using namespace std;

// -------------------------------------------------------------------------- //
// Programmer: Ryan D. Crawford
// Date: 03/25/2018
// This The sorted list class does not store any data directly. Instead,
// it contains a collection of LinkedNodeClass objects, each of which
// contains one element.
// -------------------------------------------------------------------------- //

#ifndef _SORTEDLISTCLASS_H_
#define _SORTEDLISTCLASS_H_
class SortedListClass
{
  private:
    LinkedNodeClass *head; //Points to the first node in a list, or NULL
    //if list is empty.
    LinkedNodeClass *tail; //Points to the last node in a list, or NULL
    //if list is empty.

  public:
    //Default Constructor. Will properly initialize a list to
    //be an empty list, to which values can be added.
    SortedListClass();

    //Dtor
    ~SortedListClass();

    //Clears the list to an empty state without resulting in any
    //memory leaks.
    void clear();

    //Input: The value to insert into the list
    //Allows the user to insert a value into the list. Since this
    //is a sorted list, there is no need to specify where in the list
    //to insert the element. It will insert it in the appropriate
    //location based on the value being inserted. If the node value
    //being inserted is found to be "equal to" one or more node values
    //already in the list, the newly inserted node will be placed AFTER
    //the previously inserted nodes.
    void insertValue(const int &valToInsert);

    //Prints the contents of the list from head to tail to the screen.
    //Begins with a line reading "Forward List Contents Follow:", then
    //prints one list element per line, indented two spaces, then prints
    //the line "End Of List Contents" to indicate the end of the list.
    void printForward() const;

    //Prints the contents of the list from tail to head to the screen.
    //Begins with a line reading "Backward List Contents Follow:", then
    //prints one list element per line, indented two spaces, then prints
    //the line "End Of List Contents" to indicate the end of the list.
    void printBackward() const;

    //Removes the front item from the list and returns the value that
    //was contained in it via the reference parameter. If the list
    //was empty, the function returns false to indicate failure, and
    //the contents of the reference parameter upon return is undefined.
    //If the list was not empty and the first item was successfully
    //removed, true is returned, and the reference parameter will
    //be set to the item that was removed.
    bool removeFront(int &theVal);

    //Removes the last item from the list and returns the value that
    //was contained in it via the reference parameter. If the list
    //was empty, the function returns false to indicate failure, and
    //the contents of the reference parameter upon return is undefined.
    //If the list was not empty and the last item was successfully
    //removed, true is returned, and the reference parameter will
    //be set to the item that was removed.
    bool removeLast(int &theVal);

    //Returns the number of nodes contained in the list.
    int getNumElems() const;

    //Provides the value stored in the node at index provided in the
    //0-based "index" parameter. If the index is out of range, then outVal
    //remains unchanged and false is returned. Otherwise, the function
    //returns true, and the reference parameter outVal will contain
    //a copy of the value at that location.
    int operator[]( const unsigned int index) const;
};
#endif
