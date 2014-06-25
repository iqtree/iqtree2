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
  
  /* TODO ModelSubst::rates */

  virtual void init(const char *model_name, string model_params, StateFreqType freq, string freq_params);

  void initMoranWithMutation();

  virtual int getNDim();

private:

};

#endif /* _MODELPOMO_H_ */
