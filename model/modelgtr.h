//
//modelgtr.h
//Created by James Barbetti on 06-Jul-2021
//
#include "modelmarkov.h"
#include "utils/eigendecomposition.h"

class ModelGTR : public ModelMarkov {
    public:
        typedef ModelMarkov super;
        ModelGTR(PhyloTree *tree, bool count_rates);
        ~ModelGTR() = default;

        virtual void        saveCheckpoint();
        virtual void        restoreCheckpoint();
        virtual std::string getName()       const;
        virtual std::string getNameParams() const;
        void                getNameParamsFreq(std::ostream &retname) const;
        void                init(const char *model_name, string model_params,
                                 StateFreqType freq, string freq_params,
                                 PhyloTree* report_to_tree);
        virtual void        writeInfo(ostream &out);
        virtual void        getRateMatrix(double *rate_mat);
        virtual void        setRateMatrix(double* rate_mat);
        virtual void        setStateFrequency(double* freq);
        virtual void        getQMatrix(double *q_mat);
        virtual int         getNDim() const;
        virtual int         getNDimFreq() const;
        virtual bool        scaleStateFreq();
        virtual void        setVariables(double *variables);
        virtual bool        getVariables(const double *variables);
        virtual double      targetFunk(double x[]);
        virtual bool        isUnstableParameters();
        virtual void        setBounds(double *lower_bound, double *upper_bound, 
	                                  bool *bound_check);
        virtual double      optimizeParameters(double gradient_epsilon,
                                               PhyloTree* report_to_tree);
        virtual void        decomposeRateMatrix();
        virtual void        readRates(istream &in);
        virtual void        readRates(std::string str);
        virtual void        readStateFreq(istream &in, PhyloTree* report_to_tree);
        virtual void        readStateFreq(string str, PhyloTree* report_to_tree);
        virtual void        readParameters(const char* file_name, 
                                           bool        adapt_tree_ignored,
							               PhyloTree*  report_to_tree);
        virtual void        freeMem();

        virtual void    computeTransMatrix(double time, 
                                   double *trans_matrix, int mixture=0);
        void    computeTransMatrixFreq(double time, 
                                       double* trans_matrix);
        double  computeTrans(double time, int state1, int state2);
        double  computeTrans(double time, int state1, int state2, 
                             double &derv1, double &derv2);
        virtual void computeTransDerv(double  time, double* trans_matrix, 
	                             double* trans_derv1, double* trans_derv2,
                                 int mixture = 0);
        void    computeTransDervFreq(double time, double rate_val, 
                                     double* trans_matrix, double* trans_derv1, 
				 					 double* trans_derv2);

    protected:
        bool       half_matrix;
        int        highest_freq_state;
};