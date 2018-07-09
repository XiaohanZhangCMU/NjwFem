
// The declaration vector<int[Nvar] > is illegal, 
// so use vector<intArray<Nvar> >

#ifndef __INTARRAY_H__
#define __INTARRAY_H__

template<const int Nsize>
class intArray
{
 public:
  intArray() {return;}
  
  operator int*() {return ii;}
  operator const int*() const {return ii;}  // for use in const expressions
  
 private:
  int ii[Nsize];
};
#endif                          // !defined(__INTARRAY_H__)
