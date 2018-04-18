#include <iostream>
#include <cstdlib>
#include "LinkedNodeClass.h"
#include "SortedListClass.h"
using namespace std;

// -------------------------------------------------------------------------- //
// Programmer: Ryan D. Crawford
// Date: 03/25/2018
// This The sorted list class does not store any data directly. Instead,
// it contains a collection of LinkedNodeClass objects, each of which
// contains one element.
// -------------------------------------------------------------------------- //



//Default Constructor. Will properly initialize a list to
//be an empty list, to be an empty list, to which values can be added.
SortedListClass::SortedListClass()
{
  // Create an empty list where the head and the tail are null pointers
  head = new LinkedNodeClass(0, 0, 0);
  tail = new LinkedNodeClass(0, 0, 0);
}

//Dtor
SortedListClass::~SortedListClass()
{
   clear(); //clear the list
   //Deallocate the memory allocated for head and the tail
  delete head;
  delete tail;
}

//Clears the list to an empty state without resulting in any
//memory leaks.
void SortedListClass::clear()
{
  if (getNumElems() > 0) //Make sure there is actually something to delete...
  {
    //Create a tempory pointer to store the next pointer to be deleted
    LinkedNodeClass *temp = head->getNext();

    //Until the last element is reached, update the temp variable to the
    do
    {
      //Go to the next element in the list
      temp = temp->getNext();

      //Now that the pointer has been updated to the next node, delete the
      //previous node
      LinkedNodeClass *nodeToDelete = temp->getPrev();
      delete nodeToDelete;
    }
    while(!(temp->getNext() != tail));

    //Set the head and tail pointers to null
    head->setNextPointerToNull();
    tail->setPreviousPointerToNull();
  }
}

//Input: The value to insert into the list
//Allows the user to insert a value into the list. Since this
//is a sorted list, there is no need to specify where in the list
//to insert the element. It will insert it in the appropriate
//location based on the value being inserted. If the node value
//being inserted is found to be "equal to" one or more node values
//already in the list, the newly inserted node will be placed AFTER
//the previously inserted nodes.
void SortedListClass::insertValue(const int &valToInsert)
{
  //If there are zero elements in the list, create a new LNC and make the prev
  //and next nodes the head and tail
  if (getNumElems() == 0)
  {
    //Create a new LNC object with the previous and next pointers to the
    //head and tail
    LinkedNodeClass* newNode = new LinkedNodeClass(head, valToInsert, tail);

    //Set the head and tail pointers to "this" node
    newNode->setBeforeAndAfterPointers();
  }
  else
  {
    //Create a temporay variable to store the values and pointers
    LinkedNodeClass *temp = head;

    //Wile the temp val, is greater than the next value, update the temp pointer
    //to the next node
    while(valToInsert > temp->getNext()->getValue() && temp->getNext() != tail)
    {
      temp = temp->getNext();
    }

    //Create a new "LinkedNodeClass" object
    LinkedNodeClass *newNode = new LinkedNodeClass(temp, valToInsert,
                                                   temp->getNext());

    //Set the previous and next nodes to this node
    newNode->setBeforeAndAfterPointers();
  }
}

//Prints the contents of the list from head to tail to the screen.
//Begins with a line reading "Forward List Contents Follow:", then
//prints one list element per line, indented two spaces, then prints
//the line "End Of List Contents" to indicate the end of the list.
void SortedListClass::printForward() const
{
  //Print the header for the list
  cout << "Forward List Contents Follow:" << endl;

  //Print the contents of the list if it isnt empty
  if (head->getNext() != 0)
  {
    //Create a temporary pointer to the head to the list
    LinkedNodeClass *temp = head->getNext();

    //Go to the next element in the list and print the value
    while(temp != tail)
    {
      //Output the value contained in the temp variable
      cout << "  " << temp->getValue() << endl;

      temp = temp->getNext(); // make the temporary pointer the next pointer
    }
  }
  //Mark the end of the list:
  cout << "End Of List Contents" << endl;
}

//Prints the contents of the list from tail to head to the screen.
//Begins with a line reading "Backward List Contents Follow:", then
//prints one list element per line, indented two spaces, then prints
//the line "End Of List Contents" to indicate the end of the list.
void SortedListClass::printBackward() const
{
  //Print the header for the list
  cout << "Backward List Contents Follow:" << endl;

  if (tail->getPrev() != 0)
  {
    //Create a temporary pointer to the head to the list
    LinkedNodeClass *temp = tail->getPrev();

    //Go to the next element in the list and print the value
    while(temp != head)
    {
      //Output the value contained in the temp variable
      cout << "  " << temp->getValue() << endl;

      temp = temp->getPrev(); // make the temporary pointer the next pointer
    }
  }
  //Mark the end of the list:
  cout << "End Of List Contents" << endl;
}

//Removes the front item from the list and returns the value that
//was contained in it via the reference parameter. If the list
//was empty, the function returns false to indicate failure, and
//the contents of the reference parameter upon return is undefined.
//If the list was not empty and the first item was successfully
//removed, true is returned, and the reference parameter will
//be set to the item that was removed.
bool SortedListClass::removeFront(int &theVal)
{
  if (getNumElems() == 0)
  {
    return false;
  }
  else if (getNumElems() == 1) //If there is only one element, clear the list
  {
    //Updata the value of the element being deleted
    theVal = head->getNext()->getValue();
    
    clear(); // clear the list
    return true;
  }
  else // Create a new version of the second element and delete the old first
  {    // and second elements
    //Create temporary pointer for the node to delete and the entry after
    LinkedNodeClass *fisrtNode = head->getNext();

    //Update the value to be deleted
    theVal = fisrtNode->getValue();

    LinkedNodeClass *secondNode = fisrtNode->getNext();

    //Create a new version of the second Node pointing at the third element and
    //the header
    LinkedNodeClass *newFirstNode = new LinkedNodeClass(head,
                                                        secondNode->getValue(),
                                                        secondNode->getNext());
    //Delete the old first and second nodes
    delete fisrtNode;
    delete secondNode;

    //Set the before and after pointers to this node
    newFirstNode->setBeforeAndAfterPointers();
    return true;
  }
}

//Removes the last item from the list and returns the value that
//was contained in it via the reference parameter. If the list
//was empty, the function returns false to indicate failure, and
//the contents of the reference parameter upon return is undefined.
//If the list was not empty and the last item was successfully
//removed, true is returned, and the reference parameter will
//be set to the item that was removed.
bool SortedListClass::removeLast(int &theVal)
{
  if (getNumElems() == 0)
  {
    return false;
  }
  else if (getNumElems() == 1) //If there is only one element, clear the list
  {
    //Update the value to the value about to be deleted
    theVal = tail->getPrev()->getValue();

    //Delete the node
    LinkedNodeClass *nodeToDelete = head->getNext();
    delete nodeToDelete;

    //Set the head and tail pointers to null
    head->setNextPointerToNull();
    tail->setPreviousPointerToNull();

    return true;
  }
  else // Create a new version of the second element and delete the old first
  {    // and second elements
    //Create temporary pointer for the node to delete and the entry after
    LinkedNodeClass *lastNode = tail->getPrev();

    //Update the value to be deleted
    theVal = lastNode->getValue();

    LinkedNodeClass *secondtoLastNode = lastNode->getPrev();

    //Create a new version of the last Node pointing at the third to last
    //element and the tail
    LinkedNodeClass *newLastNode = new LinkedNodeClass(
                                          secondtoLastNode->getPrev(),
                                          secondtoLastNode->getValue(),
                                          tail);

    //Delete the old first and second nodes
    delete lastNode;
    delete secondtoLastNode;

    //Set the before and after pointers to point at this node
    newLastNode->setBeforeAndAfterPointers();
    return true;
  }
}


//Returns the number of nodes contained in the list.
int SortedListClass::getNumElems() const
{
  int count = 0; //Create a counter to store the number of nodes

  if (head->getNext() != 0)
  {
    //Create a temporary pointer to the head to the list
    LinkedNodeClass *temp = head;

    //Update the count until the tail is reached
    while(temp->getNext() != tail  )
    {
      temp = temp->getNext(); // make the temporary pointer the next pointer
      count ++; // increment the counter
    }
  }
  return count;
}

//Provides the value stored in the node at index provided in the
//0-based "index" parameter. If the index is out of range, then outVal
//remains unchanged and false is returned. Otherwise, the function
//returns true, and the reference parameter outVal will contain
//a copy of the value at that location.which values can be added.
int SortedListClass::operator[]( const unsigned int index) const
{
  int outVal = -1; // value to output to the
  //If the index is invalid, return false
  if (index >= getNumElems() || index < 0)
  {
    return outVal;
  }
  else
  {
    int count = 0; //Create a counter to store the number of nodes

    //Create a temporary pointer to the head to the list
    LinkedNodeClass *temp = head->getNext();

    //Update the count and keep going to the next pointer until the index is
    //reached

    while(index != count)
    {
      temp = temp->getNext();  // make the temporary pointer the next pointer
      count ++;                // increment the counter
    }

    outVal = temp->getValue(); // update the outVal to the temps current value
    return outVal;
  }
}
