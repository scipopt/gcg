/*
 * class_instanciatated.h
 *
 *  Created on: Nov 19, 2012
 *      Author: bergner
 */

#ifndef CLASS_INSTANCIATATED_H_
#define CLASS_INSTANCIATATED_H_

template <typename T>
struct Instanciated
{
private:
   static SCIP_Bool instanciated;
public:
    Instanciated()
    {
    }

    ;
    void instanciate()
    {
   //     instanciated = TRUE;
    }

    SCIP_Bool  isInstanciated() const
    {
       return instanciated;
    };
protected:
    virtual ~Instanciated() { } ;
};

template <typename T> SCIP_Bool Instanciated<T>::instanciated( FALSE );


#endif /* CLASS_INSTANCIATATED_H_ */
