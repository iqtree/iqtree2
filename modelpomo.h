#ifndef _MODELPOMO_H_
#define _MODELPOMO_H_

#include "gtrmodel.h"

class ModelPoMo : public GTRModel
{     
  
public:  
  /** 
   * Constructor
   * 
   * @param tree associated tree for the model
   * 
   */
  ModelPoMo(PhyloTree *tree, bool count_rate = true);
  
  /* TODO ModelSubst::Rates */

  /* TODO ModelSubst::State_Frequ */

  /* TODO */
  virtual void init();

  /* TODO */
  virtual int getNDim();

  /* TODO */
  virtual void setVariables();

  /* TODO */
  virtual void getVariables();

  /** 
   * Destructor
   * 
   * 
   * @return 
   */
  ~ModelPoMo();

private:

};

#endif /* _MODELPOMO_H_ */
