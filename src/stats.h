/*
 * primitive_stats.h
 *
 * Copyright (C) 2013 Louis-Francois Handfield
 * e-mail: lfhandfield@gmail.com
 *
 * This program is free software; upon notification by email to the licensor
 * of the licencee identity and nature of use, the licencee can redistribute
 * this program and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation; either version 2
 * of the License, or (at the licencee option) any later version. As such,
 * no further notifications are required.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */


	template<class INPUT_DECL>
	class ProbabilityDistribution{
	public:
		typedef INPUT_DECL DISTRIB_INPUT_TYPE;
		static const bool IsPOD = false;
		static const bool IsComplex = false;
		static const bool NeedsAddLink = false;

		virtual double LL(const INPUT_DECL&) const=0;

		// *** critical section start ***

		// EMStep 1
		virtual void EMinit()=0;
		virtual void EMAlphainit(double)=0;

		virtual unsigned int NB_prestep() const{return 0;}

		virtual void EMpreregist(unsigned int step, const INPUT_DECL &instance, const double prob = 1.0f){}
		virtual void EMprefinit(unsigned int step){}



		void EMsaferegist(const INPUT_DECL &instance, const double prob = 1.0f);
		// EMStep 2
		virtual void EMregist(const INPUT_DECL &instance, const double prob = 1.0f)=0;
		// having the marginal probability that the observation comes from this distribution
		// record sufficient information for maximum likelihood updates, or gradient accent or heuristic updates
		// needs to also store the non-update loglikelihood

		// EMStep 3
		virtual double EMfinit()=0; // ***MAY ONLY BE CALLED ONCE, AFTER AN EMinit function!*** (exits critical section)
		// get complete likelihood, prior to the M step for gradient accent or heuristic updates updates,
		// but after the M step for distribution with maximum likelihood updates that enables us to track likelihood changes

		// *** critical section end ***


		// EMStep 4 (one of the two is to be called)

		virtual ProbabilityDistribution<INPUT_DECL>* EMnext(double alpha) const =0; // updates parameters, AND store the required information to revert the change, if needed

		//	virtual void EMfinit(ProbabilityDistribution &target, double alpha) const{} // updates parameters, AND store the required information to revert the change, if needed

		//		virtual void* EMfinit(ProbabilityDistribution &target, double alpha)=0; // construct the class instead... to get a pointer


		//	virtual double EMloglikelihood()=0; // recomputes likelihood after EM step (assumes EMfinit has been called)

		virtual void EMclear()=0;
	};


	/*
	 virutal class probability distribution2:

	 class ProbDist{



	 // EMStep 1
	 virtual void EMinit()=0;
	 virtual void EMAlphainit(double)=0;
	 void EMsaferegist(const INPUT_DECL &instance, const double prob = 1.0f);
	 // EMStep 2
	 virtual void EMregist(const INPUT_DECL &instance, const double prob = 1.0f)=0;

	 // EMStep 3
	 virtual double EMfinit(); // returns likelihood, after altering parameters which yeild a known change to the likelihood

	 // EMStep 4
	 void EMnext(ProbDist& double alpha) const; //


	 }
	 */



template<unsigned int DIMS =0u>
class DirichletDistribution{
    public:
    Tuple<double, DIMS> param;
    double konstant;
    void operator()(double& _out, Tuple<double, DIMS>& _in )const;
    double getLL(const DataGrid<double, 2u> &samples)const;
    void wrLLDerivative(const DataGrid<double, 2u> &samples, Tuple<double, DIMS> &fout)const;
    void updateKonstant();
    DirichletDistribution<DIMS>& toOne();
    void setSize(uint32_t _size){param.setSize(_size);}
    // void wrRandSample(Tuple< Tuple<double, DIMS> > & _fout);
};

class HiddenIndexArray{
public:
    Tuple<uint32_t> mainindexes;
    DataGrid<double> data;
};

template<unsigned int DIMS =0u>
class ParabolicMixingDistribution : public Oper2< double, Tuple<double, DIMS> >{
public:
    // the probability density is a simple second degree polynomial.
    // domain is constrained so "1 = \sum X_i" and "X_i >= 0 for all i", similar to DirichletDistribution
    Trianglix<double, DIMS> shape;
    Tuple<double, DIMS> minimum;
    double LL_constant;
    void setSize(int _size);
    ParabolicMixingDistribution& toOne(); // makes a uniform distribution
    void operator()(double& _out, Tuple<double, DIMS>& _in )const;
    double getLL(const DataGrid<double, 2u> &samples)const;

    template<Tuple_flag TF> double expectStep(KeyElem<uint32_t,Tuple<double, 0u, TF> > &_out, const KeyElem<uint32_t,Tuple<double, 0u, TF> > &_in, const Tuple<double, 0u, TF> &ll_deriv, const Trianglix<double, 0u> &ll_dblder, const Tuple<double> &ll_jumps) const; // update value and output LL difference

    void updateKonstant();
};



class HiddenStates{
public:
    Trianglix<double> prior;
    SparseMatrix<double> value;
    double getLogError();
};


// non-parametric and hierarchical clustering hybrid,
// hidden_state H
template<class C, unsigned int R_S =0u, unsigned int C_S =0u>
class BiClassifier{
public:
    // high-dimensional observations are made for a large number of entities, we presume that observations are explained by a pair
    // of hidden state associated to the observation and entity respectively. The objective is to maximize the likelihood of such
    // observations by Expected Maximization.

    // In addition, some additional hidden variable can be added that are specific to single entity/observation; however, these
    // should be used to normalize data only.

    // Since pairs of hidden states influence observations, they need to formally explain every combination of states.


    class Shared_EMSCOPE{
        void init_routine(BiClassifier<C,R_S,C_S>& _target);
    public:
        BiClassifier<C,R_S,C_S>& target;

        TMatrix<uint32_t> maxgrid_nonzero_count; // (par_RM.sizes[0])

        DataGrid<uint32_t> row_nbmax_hcol; // (par_RM.sizes[1] x data.size[0])
        DataGrid<uint32_t> col_nbmax_hrow; // (par_RM.sizes[0] x data.size[1])

        Tuple<uint32_t> nb_maxrowh; //(par_RM.sizes[0])
        Tuple<uint32_t> nb_maxcolh; //(par_RM.sizes[1])

        Tuple< Tuple<double ,2u> > center_shift;

        SparseMatrix<double> introns;
        SparseMatrix<double> data;
        SparseMatrix<double> trdata;
        bool istimefor_columns;
        bool tmptmp;

        bool ignore_drop;
        double hiddenalpha;
        //bool heuristictime;

        SparseMatrix<double> newhidden; // shared output
        Tuple< Tuple<double ,2u> > newscale; // shared output
        //DataGrid<uint32_t> new_nbmax_h; // (par_RM.sizes[0] x data.size[1]) or (par_RM.sizes[1] x data.size[0])
        Shared_EMSCOPE(BiClassifier<C,R_S,C_S>& trg, const SparseMatrix<C>& data,const Vector<uint32_t> &excl_list);
        Shared_EMSCOPE(BiClassifier<C,R_S,C_S>& trg, const SparseMatrix<Tuple<uint32_t, 2u> >& data, const Vector<uint32_t> &excl_list, bool incl_introns = true); // Exonic and Intronic!
   //     template<int F> void initCounts( SparseMatrix<double, 0u> &rowhid, SparseMatrix<double, 0u> &colhid);
    };
    //


    class Parallel_EMSCOPE : public Event<void>{ // const BiClassifier<C,R_S,C_S>&
    public:
        Vector<uint32_t> indexes_changes;
        TMatrix< Tuple<double,2u> > dll_dMR2;
        Trianglix< double > dd_massive;

        bool use_stdout;
        double alpha;

        double ll[2];
        Tuple<double>* rowdata;
        Tuple<double>* coldata;
        uint32_t* cachind;
        uint32_t colrange[2];
        uint32_t rowrange[2];
        uint32_t nb_fail;
        uint32_t nb_insist;

        Shared_EMSCOPE* scp;

        TMatrix<int32_t> nonzero_diff;
        Tuple<int32_t> nbmax_diff;

        void setRange(uint32_t r_start, uint32_t r_max, uint32_t c_start, uint32_t c_max);
        //void setPointers(Tuple<double>* _r, Tuple<double>* _c , uint32_t* _h);
        void setSizes(Tuple<uint32_t, 2u> &coor, int nbrow);
        uint32_t operator()();
    };

    class TaskNormalize : public Event<uint32_t>{
    public:

        TaskNormalize(uint32_t nbthread);
        void operator()(uint32_t threadID);
    };

    Trianglix<double> rowhid_prior;
    double rowhid_prior_constant;
    Trianglix<double> colhid_prior;
    double colhid_prior_constant;

    SparseMatrix<double> rowhid;
    SparseMatrix<double> colhid;
    Tuple< Tuple<double ,2u> > par_colscale;
    Tuple< Tuple<double ,2u> > par_rowscale;

    Tuple< double > par_coldropterm;
    Tuple< double > par_rowdropterm;

    MPReachDistribution row_states_prior;
    MPReachDistribution col_states_prior;
    Oper2<double, DataGrid<C, 2u> >* distrib;

    TMatrix< Tuple<double ,2u> > par_RM;
    TMatrix<double> par_D;
    TMatrix<double> par_S;

    Tuple<Tuple<double> > bg_transcripts;
    Tuple< uint32_t > bg_tr_col_index;

    //~BiClassifier();

    void setDistrib(Tuple<unsigned int, 2u> dims, Oper2<double, DataGrid<C, 2u> >* distr);
    void run2D_EM(const SparseMatrix<C>& in);
    void setRandom(int nbrows, int nbcolumns, int nbrowshidden, int nbcolumnshidden);
    void setSemiRandom(int nbrows, int nbrowshidden, Vector<uint32_t> input_col);

    ERRCODE initHidden(int nbrowshidden ,int nbcolumnshidden, int nbrows, int nbcolumns, double row_sparse_size, double col_sparse_size, bool reuseRows =false, bool reuseCols =false);

    void run2D_EM(const DataGrid<C, 2u> & in);

    void initScale(const SparseMatrix<C>& in, bool doRows= true, bool doCols =true);
    template<class P> void run2D_EM_v2(const SparseMatrix<C>& in, const Vector<uint32_t> &excl_list, P& tb, uint32_t nbstep = 200u, unsigned int nbthread = 1u);
    void run2D_EM_v3(ThreadBase &tb, const SparseMatrix<C>& in, const Vector<uint32_t> &excl_list, const TMatrix<double>& soup, const Vector<uint32_t> &soupindex, uint32_t nbstep = 200u, unsigned int nbthread = 1u);
    class Task_run2D_EM_v4;
    void run2D_EM_v4(const SparseMatrix<C>& in, const Vector<uint32_t> &excl_list, uint32_t nbstep = 200u, bool use_heavy_prior=true);

    ERRCODE run2D_EM_dropout(ThreadBase &tb, Tuple< double > &fout_rowdrop, Tuple< double > &fout_coldrop, const SparseMatrix<C>& data, uint32_t nbstep = 200u);
    template<class P> void run2D_scaleEM(const SparseMatrix<C>& in, P& tb, uint32_t nbstep = 200u, unsigned int nbthread = 1u);
    void doRecenter();
    void exportTiff(const char* folderpath)const;
    void exportZscore(const char* folderpath, const SparseMatrix<uint32_t> &input)const;


    class normalizeMK2_task;
    SparseMatrix<double> normalizeMK2(const SparseMatrix<uint32_t> &input)const;
    SparseMatrix<double> normalize(const SparseMatrix<uint32_t> &input, bool renormalize_pvalue_cdf=true, SparseMatrix<double>* opt_pvalout= NULL)const;
    SparseMatrix<double> normalize_dropout(const SparseMatrix<uint32_t> &input, SparseMatrix<double>* opt_logitpvalue= NULL)const;
    SparseMatrix<double> computeDeviation(const SparseMatrix<uint32_t> &input)const;
};

class MarkovAttractor{
public:
    Tuple<double> coefs;
    TMatrix<double> transit_and_prior; // hidden_state_prior on diagonal, transit rates on off diagonals
    void populateCoefs();
    Tuple<double> computeLL(const SparseMatrix<double>& hidden)const;
};

template<class C, unsigned int nbstates>
class Classifier : public Oper2< Tuple<double, nbstates> , C> {
public:
    static const bool IsPOD = false;
    static const bool IsComplex = false;
    static const bool NeedsAddLink = false;
    Oper2<double, C>* classes[nbstates];
    double* unknown_scope;
    double relerr;
    double antioss;
    Classifier(): unknown_scope(NULL){}
    void operator()(Tuple<double, nbstates>& _out, C& _in )const;
    void setUnknownProbability(double); // used for the Maximization updates only, the marjinal is outputed the Expectation
    double UnknownProbability(C &instance, const Tuple<double, nbstates>& prob);
    Oper2<double, C>* & operator[](const int& i) const{ return(classes[i]);}
    Oper2<double, C>* & operator[](const int& i) { return(classes[i]);}
    void EMinit();
    void EMAlphainit(double);
    void EMregist( C &instance, const Tuple<double, nbstates> prob, bool should_learn = true);
    void EMclear();
    double EMfinit();
    Classifier<C, nbstates> EMfinit(double alpha);
    void performEM(Vector<C> &observations, int nbstep);
};

	/* Classifier model
	   Allows an Outlier class, which will be not visible from outside the class,
	   its outiler frequenncy in data will be modeled my a fixed prior distribution, of the form 1 - (X / (4*mean))^2
	 */
	template<class C>
	class ClassifierV : public Oper2< double* , C> {
		Vector< ProbabilityDistribution<C>*  > classes;
		void update_to_defaut_transition();
		// routine used by HMM EMregist and ExpectHidden
		template<int dims> void directionnal_HMM_eval(int direction, const DataGrid<double, dims> &marginals, DataGrid<double, dims> &posterior, DataGrid<double, 2> &trcp) const;
		template<int dims> void directionnal_HMM_learn(int direction, const DataGrid<double, dims> &marginals, DataGrid<double, dims> &posterior, DataGrid<double, 2> &trcp);
		bool fix_unknown_dens;

		ProbabilityDistribution<C>** EMsclasses;
	public:
		static const bool IsPOD = false;
		static const bool IsComplex = false;
		static const bool NeedsAddLink = false;

		double* prior_freq;
		double* prior_freq_learn;

		double totdatabuffer;


		DataGrid<double, 2> transitions;
		bool isseq;


		unsigned int nbstates() const {return (unsigned int) classes.size();}

		void set_prior(double*, double unknown_mix = 1.0f);

		double* unknown_scope;
		double relerr;
		ClassifierV(): unknown_scope(NULL),prior_freq(NULL),isseq(false),prior_freq_learn(NULL), EMsclasses(NULL){}
		~ClassifierV(){delete[](unknown_scope);delete[](prior_freq);delete[](prior_freq_learn);}
		void operator()(double * & _out, C& _in )const;
		void operator()(double * _out, C& _in )const;



		void push_back(ProbabilityDistribution<C>*  newdist);
		ProbabilityDistribution<C>* pop_back();

		void setUnknownProbability(double); // used for the Maximization updates only, the marjinal is outputed the Expectation
		double UnknownProbability(C &instance, const double * const & prob);

		ProbabilityDistribution<C>*  & operator[](const int& i) const{ return(classes[i]);}
		ProbabilityDistribution<C>*  & operator[](const int& i) { return(classes[i]);}
		void setUnknownDensity(double);
		void setSequential(bool on_off = true);
		void setPrior(double* npr);

		void EMinit();
		void EMAlphainit(double);

		double EMLogLikelihood() const;

		unsigned int NB_prestep() const;

		void EMpregist( C &instance, const double* const prob, const double weight =1.0f);
		void EMpfinit();
		void EMregist( C &instance, const double* const prob, const double weight =1.0f);

		template<int dims> void EMregist(DataGrid< C, dims> &data, DataGrid< double, dims>* prob = NULL); //prob is a weight for each pixels, the expected hidden state are expected at the same time as the (registration) maximization step

		void EMaccept(bool doaccept = true);
		void EMnext(double alpha);

		void EMclear();
		double EMfinit();



		double performEM(Vector<C> &list, int nbstep, bool rand_init = true);

		// dimensionnal HMM!:


		template<int dims> void ExpectHidden(DataGrid< double, dims+1> &prob, const DataGrid< C, dims> &data,  bool include_unknown) const;
		// returns (P(Z | Data, param))
	};
	//


// HMM algorims and variants:

class MetaPixels{
public:
    template<class C, unsigned int DIMS> DataGrid<double,2> performEM(unsigned int nbstep, const DataGrid<C, DIMS> &obs, ClassifierV<C>&model, const DataGrid<unsigned int, DIMS> &meta, unsigned int nbmetapix, unsigned int first_val = 0);
};


template<class C, unsigned int SIZE>
	class GaussianObj : public ProbabilityDistribution< Tuple<C,SIZE> > {
	public:
        Tuple<C, SIZE> mean;
        Trianglix<C, SIZE> ihvar;
        double factor;

		GaussElem<C,SIZE> scope;


        GaussianObj(){}
        ~GaussianObj(){EMclear();}

        typename ExCo<C>::LMUL_TYPE computeProj(const Tuple<C,SIZE> &) const;
        void operator()(double &, Tuple<C,SIZE> &) const;

		void EMinit();
		void EMAlphainit(double);

		void EMregist(const Tuple<C,SIZE> &instance, const double prob= 1.0f);
		GaussianObj<C,SIZE>* EMnext(double alpha) const;
        double LL(const C&) const;

		double EMfinit();
		void EMclear();
		void show(FILE *f, unsigned int level = 0 )const;

	};

// observations 7 ~ NB(

template<unsigned int nbchannels>
class GaussianDistribution : public ProbabilityDistribution< Tuple<double, nbchannels> >{
public:
//		typedef LFHPrimitive::YESNO<true> IsPOD;
    const static bool IsPOD = false;
    static const bool NeedsAddLink = false;

    Tuple<double, nbchannels> mean;
    TMatrix<double, nbchannels, nbchannels> ihvar;
    double normfactor;
    double weight;

    double* stat_Scope;
    GaussianDistribution(): stat_Scope(NULL){}
    GaussianDistribution(const GaussianDistribution<nbchannels> &inp): mean(inp.mean), ihvar(inp.ihvar), normfactor(inp.normfactor), stat_Scope(NULL), weight(inp.weight) {}
    ~GaussianDistribution() {delete[](stat_Scope);}

//	GaussianDistribution(const GaussianDistribution<nbchannels>&);
//	const GaussianDistribution<nbchannels>& operator=(const GaussianDistribution<nbchannels>&);
    const GaussianDistribution<nbchannels>& operator=(const GaussianDistribution<nbchannels>& other);

    template<unsigned int nbout> GaussianDistribution<nbout> makeSubGaussian(const Tuple<unsigned int, nbout> & ) const;

    void setMean(const Tuple<double, nbchannels>& _mean){mean = _mean;}
    void setVariance(const Tuple<double, nbchannels>& _var);
    void setCovariance(const Trianglix<double, nbchannels>& _cov);

    void operator()(double &, Tuple<double, nbchannels> &) const;
    double LL(const Tuple<double, nbchannels> &)const;

    void EMinit();
    void EMAlphainit(double);

    void EMregist_priorcount(const GaussianDistribution &prior);

    void EMregist_priorcount_var_only(const GaussianDistribution &prior);


    void EMregist(const Tuple<double, nbchannels> &instance, const double prob = 1.0f);
    void EMclear();

    double EMfinit();
    double EMloglikelihood() const;

    ProbabilityDistribution< Tuple<double, nbchannels> >* EMnext(double alpha) const {return new GaussianDistribution<nbchannels>(*this);}

    void show(FILE* o= stdout, int level = 0) const;

    template<unsigned int out_dims> GaussianDistribution<out_dims> operator*(const TMatrix<double, nbchannels, out_dims> &) const;


    double getCorrelationCoef(unsigned int d1, unsigned int d2) const;

    TMatrix<double, nbchannels, nbchannels> getVariance(){ return ihvar.inverse() * -0.5f;}

    double getVar(int ch) {return(-0.5f / ihvar.data[ch*(nbchannels+1)]); }

    void save(FILE* f) const;
    void load(FILE* f, unsigned int size = 0);

    void setWeight(const double &w);
    void setMean_of_statistic_to_zero();
    const GaussianDistribution& operator*=(const Tuple<double,nbchannels>& scale);
};
class GaussianDistributionV : public ProbabilityDistribution< double* >{
public:
    //typedef LFHPrimitive::YESNO<true> IsPOD;
    const static bool IsPOD = false;
    unsigned int nbsizes;
    double* mean;
    DataGrid<double,2> eigenvec;
    double normfactor;
    double* stat_Scope;

    GaussianDistributionV(): nbsizes(0),stat_Scope(NULL){}
    GaussianDistributionV(const GaussianDistributionV &inp): nbsizes(inp.nbsizes),  mean(new double[inp.nbsizes]), eigenvec(inp.eigenvec), normfactor(inp.normfactor) ,stat_Scope(NULL) {}

    ~GaussianDistributionV() {delete[](stat_Scope);}

//	GaussianDistribution(const GaussianDistribution<nbchannels>&);
//	const GaussianDistribution<nbchannels>& operator=(const GaussianDistribution<nbchannels>&);
    const GaussianDistributionV& operator=(const GaussianDistributionV& other);

    void setSizes(unsigned int s);

    void setMean(const double * _mean){memcpy(mean,_mean,sizeof(double)*nbsizes);}
    void setVariance(const double * _var);

    void operator()(double &, double *  &) const;
    double LL(double* const &)const;
    void EMinit();
    void EMAlphainit(double);
    void EMregist(double *  const  & instance, const double prob = 1.0f);

    double EMfinit();
    ProbabilityDistribution< double* >* EMnext(double alpha) const {return new GaussianDistributionV(*this);}


    void EMclear();
    void show(FILE* o= stdout) const;




//		GaussianDistributionV operator+(GaussianDistributionV const& other) const; // assumes iid
//		GaussianDistributionV operator-(GaussianDistributionV const& other) const; // assumes iid
//		GaussianDistributionV operator*(Datagrid<double, 2> const& other) const; // assumes iid

    void save(FILE* f) const;
    void load(FILE* f, unsigned int size = 0);
};

class MatrixDistribution : public ProbabilityDistribution< DataGrid<double, 2u> > {
    double evalLL_routine(const SparseMatrix<double> &data, DirichletDistribution<0u>& rowprior, DirichletDistribution<0u>& colprior, DataGrid<double, 2u> &rowhidden, DataGrid<double, 2u> &colhidden, const TMatrix<double,0u,0u> &mat_M, const Trianglix<double> &mat_U, const Trianglix<double> &mat_V)const;
    double evalLL_routine(const SparseMatrix<double> &data, DataGrid<double, 2u> &rowhidden, DataGrid<double, 2u> &colhidden, const TMatrix<double,0u,0u> &mat_M, const Trianglix<double> &mat_U, const Trianglix<double> &mat_V)const;
    template<Tuple_flag TF> void expect_routine(Tuple<double, 0u, TF> &_out, const Tuple<double, 0u, TF> &_in, const double* _stat_vector, const double difsum, const Trianglix<double> innersum, const Trianglix<double> &mat_U, const uint32_t _nb, Tuple<double> &der, Tuple<double> &dir, Trianglix<double> &dder);
public:
    TMatrix<double> means;
    Trianglix<double> rowprc; // "std"
    Trianglix<double> colprc; // "std"

    double norm_constant; // {xy/2} log(2pi) - logdet(rowprc) - logdet(colprc)

    KeyElem< Tuple<Trianglix<double>, 2u>, double* >* stat_scope;

    MatrixDistribution();
    void operator()(double &dens, DataGrid<double, 2u> &input) const {dens = exp(this->LL(input));}
    void setSizes(Tuple<unsigned int, 2u> coor, bool init_with_default_prior = false);
    void getSizes(Tuple<unsigned int, 2u> &coor)const { coor[0] = rowprc.getSize(); coor[1] = colprc.getSize();}
    double LL(const DataGrid<double, 2u>&) const;

   // void setCovars(Trianglix<double>& rowcov, Trianglix<double>& colcov);


    void EMinit();
    void EMAlphainit(double){}
    void EMregist(const DataGrid<double, 2u> &instance, const double prob);

    double EMHiddenMatrix(const DataGrid<double, 2u> &data, const DataGrid<double, 2u> &rowhidden, const DataGrid<double, 2u> &colhidden);
    double EMHiddenMatrix(const SparseMatrix<double> &data, DirichletDistribution<0u>& rowprior, DirichletDistribution<0u>& colprior, DataGrid<double, 2u> &rowhidden, DataGrid<double, 2u> &colhidden);
    double EMHiddenMatrix(const SparseMatrix<double> &data, DataGrid<double, 2u> &rowhidden, DataGrid<double, 2u> &colhidden); // uses uniform prior...

    double EMfinit();
    MatrixDistribution* EMnext(double alpha) const{return new MatrixDistribution(*this);}
    void EMclear(){}
};

/*
template<class C>
class NormalDistribution{
public:
    C* scope;
    unsigned int dims;
    NormalDistribution(): dims(0){}
    NormalDistribution(unsigned int nb_dims): dims(nb_dims){scope = new C[nb_dims * nb_dims]}
    ~NormalDistribution(){if (dims) delete[](scope);}

};*/

template<unsigned int DIMS>
class EuclidianEllipse : public ProbabilityDistribution< pair< Tuple<double, DIMS> , double> >{
    bool updateGreedy;
    unsigned int fixvars;
public:
    Tuple<double, DIMS> center;
    Tuple<double, DIMS> eccent;
    double width;
    double nihvar; // neg half var inverse
    double fact; // fast coef for eval likelihood
    Tuple<double, (DIMS*4) + 3> EMscope;
    GradientSearchScope alpha_search;

    // MODES

    void setGreedyUpdates(bool input)	{updateGreedy = input;}

    void set_stddev(double val)		{nihvar = -0.5f / (val * val);fact = sqrt(nihvar / -M_PI);}
    void fix_width(bool val =true)		{fixvars = ((val) ? 1 : 0) | (fixvars &0xFFFFFFFE);}
    void fix_stddev(bool val =true)		{fixvars = ((val) ? 2 : 0) | (fixvars &0xFFFFFFFD);}
    void fix_center(bool val =true)		{fixvars = ((val) ? 4 : 0) | (fixvars &0xFFFFFFFB);}
    void fix_eccent(bool val =true)		{fixvars = ((val) ? 8 : 0) | (fixvars &0xFFFFFFF7);}

    bool has_fix_width() const 	{return (fixvars & 1);}
    bool has_fix_stddev() const 	{return (fixvars & 2);}
    bool has_fix_center() const 	{return (fixvars & 4);}
    bool has_fix_eccent() const 	{return (fixvars & 8);}


    EuclidianEllipse();
    EuclidianEllipse(const Tuple<double, 2> &in_center);
//	EuclidianEllipse(const EuclidianEllipse<DIMS>&);

    // Algerabric pivot!
    unsigned int NB_prestep() const{return(1);}

    void EMpreregist(unsigned int step, const pair< Tuple<double, DIMS> , double> &instance, const double prob = 1.0f);
    void EMprefinit(unsigned int step);

    double densityOfZscore(double Zsc) const { return fact * exp(-0.5f * Zsc*Zsc);}

    void operator()(double &, pair< Tuple<double, DIMS> , double> &) const; //

    double LL(const pair< Tuple<double, DIMS> , double> &) const;

    void EMinit();
    void EMAlphainit(double);

    void EMregist( const pair< Tuple<double, DIMS> , double> &instance, const double prob = 1.0f);
    double EMfinit();

    EuclidianEllipse<DIMS> EMnext_exec(double alpha) const; // uses the best of algebric fit
    ProbabilityDistribution< pair< Tuple<double, DIMS> , double> >* EMnext(double alpha) const{return new EuclidianEllipse<DIMS>(this->EMnext_exec(alpha));}

    double EMloglikelihood() const;
    void EMclear();
    void show(FILE* o = stdout) const;
};

	/*
	template<class A>
	class DistanceError : public ProbabilityDistribution< pair<A, double> >{
		public:
		A center;
		double radius;
		double nihvar;
		double fact;

		double inside_factor; // factor which possess inner objects

		pair<A, Tuple<double,4> >* EMscope;

		DistanceError();
		void operator()(double &, pair<A, double> &) const;

		void EMinit();
		void EMAlphainit(double);
		void EMregist( const pair<A,double> &instance, const double prob = 1.0f);

		void EMfinit();
		void* EMnext() const;
		void EMclear();
		void show(FILE* o = stdout) const;

	};
*/
class FoldedGaussianDistribution : public ProbabilityDistribution< double >{
public:
    double mean;
    double i_std; // = 1 / (sqrt2)std
    FoldedGaussianDistribution();
    FoldedGaussianDistribution(double mean, double std);
    void operator()(double &_out_prob, double &obs) const;
    double LL(const double&)const;
    void EMinit();
    void EMAlphainit(double);
    void EMregist( const double &instance, const double prob = 1.0f);
    double EMfinit();
    FoldedGaussianDistribution EMnext_exec(double alpha) const;
    ProbabilityDistribution< double >* EMnext(double alpha) const{ return new FoldedGaussianDistribution( this->EMnext_exec(alpha)  ); }
    void EMclear();
    void show(FILE* o = stdout) const;
};
class VarGammaDistribution: public ProbabilityDistribution< double >{
public:
    double lambda;
    double correlation;
    void operator()(double &_out_prob, double &obs) const;
    double LL(const double&)const;
};

	// distribution which gives equal probability density to any instance, but explains a fixed portion of the data


	template<class A, class B>
	class DistribPair : public ProbabilityDistribution< pair<typename A::DISTRIB_INPUT_TYPE, typename B::DISTRIB_INPUT_TYPE > > {
	public:
		A distA;
		B distB;
		DistribPair(){}
		DistribPair(const A &f_a, const B &f_b): distA(f_a), distB(f_b){}


		void operator()(double &_out_prob, pair<typename A::DISTRIB_INPUT_TYPE, typename B::DISTRIB_INPUT_TYPE > &obs) const{double tmp; distA(tmp,obs.first); distB(_out_prob,obs.second); _out_prob *=tmp;}
		double LL(pair<typename A::DISTRIB_INPUT_TYPE, typename B::DISTRIB_INPUT_TYPE > &obs) const{return distA.LL(obs.first) + distB.LL(obs.second); }

		void EMinit(){distA.EMinit();distB.EMinit();}
		void EMAlphainit(double f){distA.EMAlphainit(f);distB.EMAlphainit(f);}
		double EMfinit(){return distA.EMfinit() + distB.EMfinit();}
		ProbabilityDistribution< pair<typename A::DISTRIB_INPUT_TYPE, typename B::DISTRIB_INPUT_TYPE > >* EMnext(double alpha){return new DistribPair(A(*((const A*)distA.EMnext(alpha))), B(*(const B*)distB.EMnext(alpha)));}
		void EMclear(){distA.EMclear();distB.EMclear();}
		void show(FILE* o = stdout) const{distA.show(o);fprintf(o,"\n");distB.show(o);}

		void EMregist(const pair<typename A::DISTRIB_INPUT_TYPE, typename B::DISTRIB_INPUT_TYPE > &instance, const double prob = 1.0f){distA.EMregist(instance.first);distB.EMregist(instance.second);}
	};



	template<class INPUT>
	class UnknownDistribution : public ProbabilityDistribution< INPUT >{
	public:
		double frequency;
		double density;
		double* statbuf; // 4 entries,

		UnknownDistribution(const double _freq) : frequency(_freq), density(pow(2.0f, -1000.0f)),statbuf(NULL) {}
		void operator()(double &, INPUT &) const;
		void EMinit();
		void EMAlphainit(double);
		double EMLogLikelihood() const;
		void EMregist(const INPUT &instance,  const double prob = 1.0f);
		void EMclear();
		double EMfinit();

		void InitialGuess(double totalLikelyhood);
	};

	// metric should have the following function:
	// double (PARAM& ,PARAM&);


	template<unsigned int NBPARAM, unsigned int NBINPUT>
	class GaussianProcessProbSimple_Heap : public ConstGrid<double, 2> {
		public:
        double noise;
		double *scale;
		Vector< KeyElem< Tuple<double, NBPARAM> , pair< Tuple<double, NBINPUT> , double > > > registered;
		virtual double operator()(const Tuple<unsigned int, 2> &coor) const{
			if (coor[0] == coor[1]) return noise * registered[coor[1]].d.second * registered[coor[0]].d.second;
			unsigned int j,k,l;
			Tuple<double, NBPARAM> dif;
			double dist, sum;
			dif = registered[coor[1]].k - registered[coor[0]].k;
			// multiply my upper triangular and computes self-innerproduct:
			k=0;
			sum =0;
			for(j=0;j<NBPARAM;j++)	{
				dist = dif[j] * scale[k];
				for(l=j+1;l<NBPARAM;l++,k++) dist += dif[l] * scale[k];
				sum += dist*dist;
			}
			return exp(sum) * registered[coor[1]].d.second * registered[coor[0]].d.second ; // covariance !
		}
		virtual void getDims(Tuple<unsigned int, 2> &o_dims) const{ o_dims[0] = registered.size(); o_dims[0] = o_dims[1];}
	};



	template<unsigned int NBPARAM, unsigned int NBINPUT>
	class GaussianProcessProbSimple : public Oper2<double, KeyElem< Tuple<double, NBPARAM> , Tuple<double, NBINPUT> > >{
     public:
     double setDerivative_routine(double* deriv,bool showamt = false );



        // Pred = K(x,x') * K(x,x)^-1 * OBS
        // Pred = K(x,x') * K(x,x)^-1 * OBS

		//  ~ N  Param \cdot Weight | noise), where Weight is a vector with prior distribution

//		double scale[TEMPLATE_TRIANGLE_NUMBER<2, NBPARAM >::ans +1]; // lower triangular root matrix, and a signal level variable n*(n+1)/2 + 1

		CurvatureSearchScope alpha_search;
		GaussianDistribution<NBINPUT> noise;

        double param_prior_mean[NBPARAM*2+1];
        double param_prior_std[NBPARAM*2+1];

        double scope[NBPARAM*2+1];
        double deriv[NBPARAM*2+1];
        double llikelyhood;
		int EMstep_count;

        unsigned int nb_param;
        Tuple< KeyElem< Tuple<double, NBPARAM> , Tuple<double, NBINPUT> >, 0u> inner_param;

        bool resample_to_fixgrid;

		// GaussianProcessProbSimple_Heap< NBPARAM, NBINPUT> EM_instances;

        Vector< KeyElem< double, KeyElem< Tuple<double, NBPARAM> , Tuple<double, NBINPUT> > > > EM_obs;

		GaussianProcessProbSimple(): resample_to_fixgrid(false), nb_param(0){alpha_search.init(1.0f, TEMPLATE_TRIANGLE_NUMBER<2, NBPARAM >::ans +1);}
		~GaussianProcessProbSimple(){}


		void EMPreinit(double* prior_mean, double* prior_std){
				alpha_search.init(0.0000001f,NBPARAM*2+1);EMstep_count=0;memcpy(param_prior_mean,prior_mean,sizeof(double)*(NBPARAM*2+1));memcpy(param_prior_std,prior_std,sizeof(double)*(NBPARAM*2+1));
				for(unsigned int i=0;i<NBPARAM*2+1;i++) scope[i] = param_prior_mean[i]  + (param_prior_std[i] * sampleGaussian() );
			}
		void EMinit();
		void EMAlphainit(double);
		void EMregist(const KeyElem< Tuple<double, NBPARAM> , Tuple<double, NBINPUT> > &instance,  const double prob = 1.0f);
		double EMfinit();
        GaussianProcessProbSimple<NBPARAM, NBINPUT>* EMnext(double alpha) const{ GaussianProcessProbSimple<NBPARAM, NBINPUT>* fout = new GaussianProcessProbSimple<NBPARAM, NBINPUT>(); EMnext(alpha, *fout); return fout;}
        void EMnext(double alpha, GaussianProcessProbSimple<NBPARAM, NBINPUT> &fout) const;
        void EMstep(double alpha);

		void EMclear();
		double EMLogLikelihood() const{return llikelyhood;}


		void setGrid(const Vector< Tuple<double, NBPARAM> > &ingr);



		Tuple<double, NBINPUT> getSignalMean(const Tuple<double, NBPARAM>&) const; // DONE!

		Tuple<double, NBINPUT> getMean(const Tuple<double, NBPARAM> &par) const {return( getSignalMean(par) + noise.mean);}
        GaussianDistribution<NBINPUT> marginal(const Tuple<double, NBPARAM> &par) const{ GaussianDistribution<NBINPUT> fout(noise); fout.mean += getSignalMean(par); return(fout);}
		void operator()(double &fout, KeyElem< Tuple<double, NBPARAM> , Tuple<double, NBINPUT> > &data) const{fout = exp(noise.LL(data.d - getSignalMean(data.k) ));}
		double LL(KeyElem< Tuple<double, NBPARAM> , Tuple<double, NBINPUT> > &data) const{return noise.LL(data.d - getSignalMean(data.k));}

        void show(FILE* f = stdout, int level =0)const;

	};



/*
class MatrixAssociation : TMatrix<double>{
public:
    typedef INPUT_DECL FUNCTION_INPUT_TYPE;
    typedef TMatrix<double> EMSCOPE_TYPE;

    double operator()(const KeyElem<Tuple<uint32_t, 2u>, Tuple<Tuple<double, 0u, TF>, 2u>  > &argument)const;
};

template<class F, class A = F::FUNCTION_INPUT_TYPE, class G = F>
class InflatedNegativeBinomialFunction : {
public:
    F mat_r;
    F mat_p;
    G mat_d;
    class EMscope{
        typename F::EMSCOPE_TYPE mat_r_scope;
        typename F::EMSCOPE_TYPE mat_p_scope;
        typename F::EMSCOPE_TYPE mat_d_scope;
    };

    template<Tuple_flag TF> double operator()(uint32_t value, Tuple<Tuple<double, 0u, TF>, 2u>)const;

static double dLL_dr(double x, double r, double p);
static double dLL_dp(double x, double r, double p);

};*/

#undef LFHTEMP
#define LFHTEMP template<unsigned int nbchannels>


LFHTEMP const GaussianDistribution<nbchannels>& GaussianDistribution<nbchannels>::operator=(const GaussianDistribution& other){
		mean = other.mean;
		ihvar = other.ihvar;
		normfactor = other.normfactor;
		return(*this);
	}

LFHTEMP template<unsigned int nbout> GaussianDistribution<nbout> GaussianDistribution<nbchannels>::makeSubGaussian(const Tuple<unsigned int, nbout> & channels) const{GaussianDistribution<nbout> fout;
	unsigned int i,j;



	TMatrix<double, nbchannels, nbchannels> varmat = ihvar.inverse();
	TMatrix<double, nbout, nbout> var_new;


	for(i=0;i<nbout;i++){
		fout.mean[i] = mean[channels[i]];
		for(j=0;j<i;j++) var_new.data[i + j * nbout] = var_new.data[j + i * nbout];
		for(j=i;j<nbout;j++) var_new.data[i + j * nbout] = varmat.data[channels[i] + channels[j] * nbchannels];
	}


	DataGrid<double,2> dainv;
	Vector<double> eigen;

	dainv.makeDiagonalizer_ofinverse(var_new,eigen);
	TMatrix<double, nbout, nbout > funt(dainv);

	var_new = funt.mkTranspose();

	for(i=0;i<nbout;i++) for(j=0;j<nbout;j++) funt.data[i+j*nbout] *= eigen[i];


	funt *= var_new;

	fout.normfactor = eigen[0];for(i=1;i<nbout;i++) fout.normfactor *= eigen[i];
	fout.normfactor = sqrt(fabs(fout.normfactor) / M_PI);

	fout.ihvar = funt;


	return(fout);
	}

LFHTEMP template<unsigned int out_dims> GaussianDistribution<out_dims> GaussianDistribution<nbchannels>::operator*(const TMatrix<double,nbchannels, out_dims> &mat) const{

	// mean <- mat * mean
	// invar <- mat * invar * mat

	GaussianDistribution<out_dims> fout;

	fout.mean = mat * mean;
	fout.weight = weight;

	// TODO


	return(*fout);
}


LFHTEMP void GaussianDistribution<nbchannels>::EMinit(){
		if (stat_Scope == NULL) stat_Scope = new double[2 + ((nbchannels * (nbchannels +3)) / 2) ];
		memset(stat_Scope, '\0', sizeof(double)*  (2 + ((nbchannels * (nbchannels +3)) / 2)) );

	}

LFHTEMP void GaussianDistribution<nbchannels>::EMAlphainit(double alpha){
		if (stat_Scope == NULL) {stat_Scope = new double[2 + ((nbchannels * (nbchannels +3)) / 2) ];
		memset(stat_Scope, '\0', sizeof(double)*  (2 + ((nbchannels * (nbchannels +3)) / 2)) );
		}else{
			stat_Scope[0] *= (1.0f - alpha);
			stat_Scope[1] *= (1.0f - 2.0f * alpha + alpha*alpha);
			int i,k; k = 2 + ((nbchannels * (nbchannels +3))/ 2);
			for(i=2;i<k;i++) stat_Scope[i] *= (1.0f - alpha);
		}
	}

LFHTEMP void GaussianDistribution<nbchannels>::EMregist( const Tuple<double, nbchannels> &instance, const double prob){
		if ((prob <= 0.0f)||(!ExOp::isValid(prob))) return;
		int i;
		for(i=0;i<nbchannels;i++) if (!ExCo<double>::isValid(instance[i])) break;

		if (i != nbchannels){
			printf("warning! NAN in input data for gaussian!\n");
			return;
		}
		stat_Scope[0] += prob;
		stat_Scope[1] += prob * prob;
		int j,k;
		k=2;

		for(i=0;i<nbchannels;i++,k++) stat_Scope[k] += instance[i] * prob ;
		for(i=0;i<nbchannels;i++){
			for(j=i;j<nbchannels;j++, k++) stat_Scope[k] += instance[i] * instance[j] * prob ;
		}


	//	LFHDebug_RETURN_void(((nbchannels == 1)&&(stat_Scope[0] * stat_Scope[3] - stat_Scope[2] * stat_Scope[2] < 0.0f)));

	//	printf("warning! neg variance! %e\n", stat_Scope[0] * stat_Scope[3] - stat_Scope[2] * stat_Scope[2]);
	//	printf("%e\t%e\t%e\t%e\n", stat_Scope[0],  stat_Scope[1], stat_Scope[2], stat_Scope[3]);



	}

LFHTEMP void GaussianDistribution<nbchannels>::EMregist_priorcount(const GaussianDistribution<nbchannels> &instance){
	for(unsigned int i= 1 + (((nbchannels + 3) * nbchannels) / 2);i != 0xFFFFFFFF;i--) stat_Scope[i] += instance.stat_Scope[i];
}

LFHTEMP void GaussianDistribution<nbchannels>::EMregist_priorcount_var_only(const GaussianDistribution<nbchannels> &instance){
	unsigned int i,j,k;
	if (instance.stat_Scope == NULL) exit(123);
	if (stat_Scope == NULL) exit(123);

	if (stat_Scope[0] <  ExCo<double>::epsilon()){
		for(i= 1 + (((nbchannels + 3) * nbchannels) / 2);i != 0xFFFFFFFF;i--) stat_Scope[i] = instance.stat_Scope[i];
	}else{
	double factor = instance.stat_Scope[0] /  (stat_Scope[0] * stat_Scope[0]);
	k = nbchannels+2;
	for(i=0;i<nbchannels;i++){
		for(j=i;j<nbchannels;j++, k++) stat_Scope[k] += instance.stat_Scope[k] + factor * stat_Scope[2+i] * stat_Scope[2+j] ;
	}

	factor = 1.0f + (instance.stat_Scope[0] / stat_Scope[0]);

	for(k=2;k<nbchannels+2;k++) stat_Scope[k] *= factor;

	stat_Scope[0] += instance.stat_Scope[0];
	stat_Scope[1] += instance.stat_Scope[1];
	}
}

LFHTEMP void GaussianDistribution<nbchannels>::setWeight(const double &w){
		double factor = w / stat_Scope[0];
		stat_Scope[0] = w;
		stat_Scope[1] *= factor * factor;
		for(unsigned int i= 1 + (((nbchannels + 3) * nbchannels) / 2);i > 1 ;i--) stat_Scope[i] *= factor;
	}

LFHTEMP void GaussianDistribution<nbchannels>::setMean_of_statistic_to_zero(){
	unsigned int i,j,k;
	for(k=2;k<nbchannels+2;k++) stat_Scope[k] /= sqrt(stat_Scope[0]);

	for(i=0;i<nbchannels;i++){
		for(j=i;j<nbchannels;j++,k++) stat_Scope[k] -= stat_Scope[i+2] * stat_Scope[j+2];
	}

	for(k=2;k<nbchannels+2;k++) stat_Scope[k] = 0.0f;
}
/*

	 LFHTEMP double GaussianDistribution<nbchannels>::EMfinit(){
	 int i,j,k;
	 TMatrix<double, nbchannels, nbchannels> nihvar_old;

     Trianglix<double, 0u> nihvar; nihvar.setSize(nbchannels);
	 weight = stat_Scope[0];
	 //double num = (stat_Scope[1] - stat_Scope[0]*stat_Scope[0]) /2;
	 double num = -0.5f * stat_Scope[0]*stat_Scope[0];
	 double tmp = 1.0f / stat_Scope[0];
	 for(i=0;i<nbchannels;i++) {
	 mean[i] = stat_Scope[i+2] * tmp ;
	 if (!ExCo<double>::isValid(mean[i])) ExOp::toZero(mean[i]);
	 }
		 ExOp::toZero(nihvar);

	 k = 2 + nbchannels;
	 for(i=0;i<nbchannels;i++){
	 tmp = stat_Scope[0] * stat_Scope[k] - stat_Scope[2+i] * stat_Scope[2+i];
	 if (tmp <= 0) tmp = DBL_EPSILON * (fabs(stat_Scope[0] * stat_Scope[k]) +  fabs(stat_Scope[2+i] * stat_Scope[2+i]));
	 nihvar.cell(i,i) = tmp;
	 nihvar_old.data[i * (nbchannels+1)] = tmp;
	 k++;
	 for(j=i+1;j<nbchannels;j++, k++){
	 tmp = stat_Scope[0] * stat_Scope[k] - stat_Scope[2+i] * stat_Scope[2+j];
	 if (tmp <= 0) tmp = DBL_EPSILON * (fabs(stat_Scope[0] * stat_Scope[k]) +  fabs(stat_Scope[2+i] * stat_Scope[2+j]));
	 nihvar.cell(j,i) = tmp;
		 nihvar_old.data[j + i * nbchannels] = nihvar_old.data[i + j * nbchannels] = tmp;
	 }
	 }

		 nihvar_old /= num;

		 //nihvar.show();

		 DataGrid<double,2> dainv;
		 Vector<double> eigen;

		 dainv.makeDiagonalizer_ofinverse(nihvar_old,eigen);


		 TMatrix<double, nbchannels, nbchannels> funt(dainv);



		 nihvar_old = funt.transpose();

		 for(i=0;i<nbchannels;i++) for(j=0;j<nbchannels;j++) funt.data[i+j*nbchannels] *= eigen[i];


		 funt *= nihvar_old;

		 Trianglix<double, 0u> nihvar2 = nihvar.inverse(); // _MK3
		// nihvar2.show();
		// funt.show();
		 double tmp_old;
		 tmp_old = eigen[0];for(i=1;i<nbchannels;i++) tmp_old *= eigen[i];
		 tmp_old = sqrt(fabs(tmp_old)) * pow(M_PI, -0.5f * nbchannels);

	// if (!ExOp::isValid(nihvar2.data[0]))
	//	 nihvar2 = nihvar.inverse();

	//	 ExOp::show( Matrix<double>(nihvar) *  Matrix<double>(nihvar2) );
	 tmp = nihvar.log_determinant();
		 printf("%e\t%e\n", tmp_old, exp( -0.5f * tmp) * pow(-num / M_PI , 0.5f * nbchannels )); //
	//	 printf("hohoho! %e\t%e\n", nihvar.log_determinant() , nihvar2.log_determinant());
		 if ( (ExOp::isValid(nihvar2)) && (ExOp::isValid(tmp)) ){


			 normfactor = exp( -0.5f * tmp) * pow(-num / M_PI , 0.5f * nbchannels );



	 for(i=0;i<nbchannels;i++){
		 ihvar.data[i * (nbchannels+1)] = nihvar2.cell(i,i) * num;
		 for(j=i+1;j<nbchannels;j++){
			 ihvar.data[j + i * nbchannels] = ihvar.data[i + j * nbchannels] = nihvar2.cell(j,i) * num;
		 }
	 }

     } else{
		 printf("NOSE!\n");
	 }



	 return( stat_Scope[0] * (log(normfactor) - 0.5f * nbchannels));

	 }
	*/

LFHTEMP double GaussianDistribution<nbchannels>::EMfinit(){
    int i,j,k;
    TMatrix<double, nbchannels, nbchannels> nihvar;
    weight = stat_Scope[0];
	double num = (stat_Scope[1] - stat_Scope[0]*stat_Scope[0]) /2;
	//	double num = -0.5f * stat_Scope[0]*stat_Scope[0];
		double tmp = 1.0f / stat_Scope[0];
		for(i=0;i<nbchannels;i++) {
		mean[i] = stat_Scope[i+2] * tmp ;
			if (!ExCo<double>::isValid(mean[i])) ExOp::toZero(mean[i]);
		}


		k = 2 + nbchannels;
		for(i=0;i<nbchannels;i++){
			tmp = stat_Scope[0] * stat_Scope[k] - stat_Scope[2+i] * stat_Scope[2+i];
			if (tmp <= 0) tmp = DBL_EPSILON * (fabs(stat_Scope[0] * stat_Scope[k]) +  fabs(stat_Scope[2+i] * stat_Scope[2+i]));
			nihvar.data[i * (nbchannels+1)] = tmp;
			k++;
			for(j=i+1;j<nbchannels;j++, k++){
				tmp = stat_Scope[0] * stat_Scope[k] - stat_Scope[2+i] * stat_Scope[2+j];
				if (tmp <= 0) tmp = DBL_EPSILON * (fabs(stat_Scope[0] * stat_Scope[k]) +  fabs(stat_Scope[2+i] * stat_Scope[2+j]));
				nihvar.data[j + i * nbchannels] = nihvar.data[i + j * nbchannels] = tmp;
			}
		}



		nihvar /= num;

		//nihvar.show();

		DataGrid<double,2> dainv;
		Vector<double> eigen;

		dainv.makeDiagonalizer_ofinverse(nihvar,eigen);


		TMatrix<double, nbchannels, nbchannels> funt(dainv);



		nihvar = funt.mkTranspose();

		for(i=0;i<nbchannels;i++) for(j=0;j<nbchannels;j++) funt.data[i+j*nbchannels] *= eigen[i];


		funt *= nihvar;




		tmp = eigen[0];for(i=1;i<nbchannels;i++) tmp *= eigen[i];
		tmp = sqrt(fabs(tmp)) * pow(M_PI, -0.5f * nbchannels);

	if ((ExOp::isValid(funt)) && (ExOp::isValid(tmp)) ) {ihvar = funt; normfactor = tmp;}
		return( stat_Scope[0] * (log(normfactor) - 0.5f * nbchannels));

	}

LFHTEMP double GaussianDistribution<nbchannels>::EMloglikelihood() const{return stat_Scope[0] * (log(normfactor) - 0.5f* nbchannels);}


LFHTEMP void GaussianDistribution<nbchannels>::EMclear(){delete[](stat_Scope); stat_Scope = NULL; }

LFHTEMP void GaussianDistribution<nbchannels>::operator()(double &_out, Tuple<double, nbchannels> &in) const{
		Tuple<double , nbchannels> diff = in - mean;
		Tuple<double , nbchannels> xform = ihvar * diff;

	//	printf("%f\t%f\t%f\n",diff[0], ihvar.data[0], xform[0]);

		double sum = diff.dotProduct(xform);



		_out = normfactor * exp(sum);

	/*	LFHDebug_RETURN_void((!ExCo<double>::isValid(_out)) || (_out < 0.0f));
		in.show(stdout);
		ihvar.show(stdout);
		printf("%e\t%e\t%e\n",_out, sum,normfactor);

		while(true) ExCo<double>::isValid(_out);
		exit(1);*/
	}

LFHTEMP double GaussianDistribution<nbchannels>::LL(const Tuple<double, nbchannels> &in) const{
    Tuple<double , nbchannels> diff = in - mean;
    Tuple<double , nbchannels> xform = ihvar * diff;

    //	printf("%f\t%f\t%f\n",diff[0], ihvar.data[0], xform[0]);

    double sum = diff.dotProduct(xform);
    return( log(normfactor) + sum );
}
LFHTEMP void GaussianDistribution<nbchannels>::setVariance(const Tuple<double, nbchannels>& _var){
    int i,j;
    normfactor = 1.0f;
    for(i=0;i<nbchannels;i++){
        for(j=0;j<nbchannels;j++){
            ihvar.data[i + j * nbchannels] = ( i != j) ? 0.0f : -2.0f * _var[i];
        }
        normfactor *= ihvar.data[i*(nbchannels+1)];
    }
    Trianglix<double> mt = Trianglix<double>(ihvar);
    normfactor = sqrt(1.0f / (M_PI*fabs(mt.determinant())));
    ihvar = TMatrix<double,nbchannels,nbchannels>(mt.mkInverse());
}
LFHTEMP void GaussianDistribution<nbchannels>::setCovariance(const Trianglix<double, nbchannels>& _cov){
    uint32_t i,j,k;
    normfactor = 1.0f;

    Trianglix<double, nbchannels> dainv = _cov.inverse();

    for(i=0;i<nbchannels;i++){
        for(j=0;j<i;j++,k++){
            ihvar.data[i + j * nbchannels] = -dainv.data[k] * 0.5;
            ihvar.data[j + i * nbchannels] = ihvar.data[i + j * nbchannels];
        }
        ihvar.data[i + j * nbchannels] = -0.5 * dainv.data[k++];
    }
    normfactor = sqrt(1.0f / (M_PI*fabs(_cov.determinant())));
    printf("Set Blur window, norm factor %f\n", normfactor);
    ihvar.show();
}
LFHTEMP double GaussianDistribution<nbchannels>::getCorrelationCoef(unsigned int d1, unsigned int d2) const{
		TMatrix<double, nbchannels, nbchannels> varmat = ihvar.inverse();
		varmat *= -0.5f;
		return( varmat.data[d1+ d2*nbchannels] / sqrt(varmat.data[d1*(1+nbchannels)]*varmat.data[d2*(1+nbchannels)] ) );
	}

LFHTEMP const GaussianDistribution<nbchannels>&  GaussianDistribution<nbchannels>::operator*=(const Tuple<double,nbchannels>  &scale){
		unsigned int i,j;

		for(i=0;i<nbchannels;i++) {
			mean[i] *= scale[i];
			normfactor /= scale[i];
			for(j=0;j<nbchannels;j++) ihvar.data[i + j*nbchannels] /= scale[i] * scale[j];
		}
		return(*this);
	}

LFHTEMP void GaussianDistribution<nbchannels>::show(FILE* o, int level) const{
		printf("nbpt=%f\tmean=", weight); ExOp::show(mean,stdout,0);
		printf("norm_factor=%e\tcovariance(upper) and correlation(lower) triangles=\n",normfactor);
		TMatrix<double, nbchannels, nbchannels> varmat = ihvar.inverse();
		varmat *= -0.5f;
		for(unsigned int i=0; i< nbchannels;i++) for(unsigned int j=i+1; j< nbchannels;j++){
			varmat.data[i+ j * nbchannels] *= pow(varmat.data[i * (1+nbchannels)] * varmat.data[j * (1+nbchannels)] ,-0.5f);
		}
		varmat.show(o);
	}

LFHTEMP void GaussianDistribution<nbchannels>::save(FILE* f) const{
		fwrite(&weight,sizeof(double),1,f);
		fwrite(&normfactor,sizeof(double),1,f);
		mean.save(f);
		fwrite(&ihvar,sizeof(TMatrix<double, nbchannels, nbchannels>),1,f);
	}
LFHTEMP void GaussianDistribution<nbchannels>::load(FILE* f, unsigned int size){
		fread(&weight,sizeof(double),1,f);
		fread(&normfactor,sizeof(double),1,f);
		mean.load(f);
		fread(&ihvar,sizeof(TMatrix<double, nbchannels, nbchannels>),1,f);
	}

	template<class INPUT>
	void UnknownDistribution<INPUT>::operator()(double & _out, INPUT &) const{
		_out = density;
	}

	template<class INPUT>
	void UnknownDistribution<INPUT>::EMinit(){
		if (statbuf == NULL) statbuf = new double[5];
		memset(statbuf,'\0',sizeof(double)*5);
	}
	template<class INPUT>
	void UnknownDistribution<INPUT>::EMAlphainit(double val){ // immune to alpha...
		if (statbuf == NULL) statbuf = new double[5];
		memset(statbuf,'\0',sizeof(double)*5);
	}
	template<class INPUT>
	void UnknownDistribution<INPUT>::EMregist(const INPUT &instance,  const double prob){
		/*statbuf[0] += 1.0f;
		statbuf[1] += prob;
		double itot = prob / density;
		if (ExCo<double>::isValid(itot)){ // the total density is too huge compared to the unknown densit
			double tmp = (1.0f - prob) * itot;
			statbuf[2] += tmp;
			statbuf[3] -= tmp* itot;
			statbuf[4] += tmp* itot* itot;
		}*/

		// linear space updates

		statbuf[0] += 1.0f;
		 statbuf[1] += prob;
		 double tmp = (1.0f - prob) * prob;
		 statbuf[2] += tmp;
		 statbuf[3] -= tmp* prob;
		 statbuf[4] += tmp* prob* prob;

		//log space updates
	/*	statbuf[0] += 1.0f;
		statbuf[1] += prob;
		double tmp = (1.0f - prob) * prob;
		statbuf[2] += tmp;
		statbuf[3] += tmp * (1 - 2 * prob);
		statbuf[4] += tmp * (6 * (prob-1) * prob + 1) ;
	 */
	}
	template<class INPUT>
	void UnknownDistribution<INPUT>::EMclear(){
		if (statbuf != NULL) {delete[](statbuf);	statbuf = NULL;}
	}

	template<class INPUT>
	double UnknownDistribution<INPUT>::EMLogLikelihood() const{
		return(0.0f);
	}

	template<class INPUT>
	double UnknownDistribution<INPUT>::EMfinit(){
		if (statbuf[0] == 0.0f) return 0.0;

		double poly[4];


		double shift;

	if ((ExCo<double>::isValid(statbuf[3]))&&(ExCo<double>::isValid(statbuf[4]))&&(statbuf[3] != 0.0f)&&(statbuf[4] != 0.0f)){
		poly[0] = statbuf[1] - frequency * statbuf[0];
		poly[1] = statbuf[2];
		poly[2] = statbuf[3];
		poly[3] = statbuf[4];
//		poly[2] = 2 * statbuf[3] + statbuf[2];
//		poly[3] = 6 * statbuf[4] + 6 * statbuf[3] + statbuf[2];
		printf("%e\t%e\t%e\t%e\n",poly[0], poly[1], poly[2] , poly[3]  );

		printf("%e true freq!\n",statbuf[1] / statbuf[0]);

		shift = -CubicRealRoot(poly,false);


		if (!ExOp::isValid(shift)){

			shift = (statbuf[1] - frequency * statbuf[0]) / statbuf[2];
			if (!ExOp::isValid(shift)) shift = sqrt(density);

		}

		}else{

			shift = (statbuf[1] - frequency * statbuf[0]) / statbuf[2];
			if (!ExOp::isValid(shift)) shift = sqrt(density);

		}

		// linear update:
		//if (shift >= 1.0f) density *= pow(2.0f,-16.0f);
		//else density *= (1.0f - shift);
		// log update:
		/*
		printf("%e\t%e\n", shift, exp(-shift) );
		if (fabs(shift) > 300.0f) shift = shift < 0.0f ? -300.0f : 300.0f;
		poly[0] = density * exp(-shift);
		while ((poly[0] == 0)||(fabs(poly[0]) == 1.0f / 0.0f)){
			printf("LOOP!!!\n");
			shift /=2;
			poly[0] = density * exp(-shift);
		}
		density = poly[0];
		*/


		if (shift < 0.0f) density *= 1.0f - shift;
		else density *= exp(-shift);
		return( density);
	}

	template<class INPUT>
	void UnknownDistribution<INPUT>::InitialGuess(double totalLikelyhood){
		density = totalLikelyhood;

	}







#undef LFHTEMP
#define LFHTEMP template<class C>


LFHTEMP	void ClassifierV<C>::setUnknownProbability(double val){
		if ((val == 0.0f)^(unknown_scope == NULL)){
			if (val == 0) {delete[](unknown_scope); unknown_scope = NULL;}
			else{
				unknown_scope = new double[10];
				unknown_scope[0] = 0.0f;
				unknown_scope[1] = val;
				unknown_scope[2] = 0.0f; // step_size
			}
		}
		fix_unknown_dens=false;
	}



LFHTEMP	void ClassifierV<C>::setUnknownDensity(double indens){
		if (unknown_scope == NULL){	unknown_scope = new double[8];}
		unknown_scope[0] = indens;
		unknown_scope[1] = 1.0f; // prior on unkwown
		fix_unknown_dens=true;
	}

LFHTEMP	double ClassifierV<C>::UnknownProbability(C &instance, const double * const & prob){
		int i;
		double sum =0.0f;
		double tmp;
		int nbstates = classes.size();
		for(i=0;i<nbstates;i++){(*classes[i])(tmp,instance);
			if (tmp != 0.0f) sum += prob[i] * log(tmp);
		}
		return(unknown_scope[0] / (exp(sum) + unknown_scope[0]));
	}



	LFHTEMP	void ClassifierV<C>::operator()(double * _out, C& _in )const{
		int i;
		double sum = 0.0f;
		int nbstates = classes.size();
		for(i=0;i<nbstates;i++){
			(*classes[i])(_out[i],_in);
			sum += _out[i];
		}

		for(i=0;i<nbstates;i++) _out[i] /= sum;

	}
LFHTEMP	void ClassifierV<C>::operator()(double * &_out, C& _in )const{
		int i;
		double sum = 0.0f;
		int nbstates = classes.size();
		for(i=0;i<nbstates;i++){
			(*classes[i])(_out[i],_in);
			_out[i] *= prior_freq[i];
			sum += _out[i];
		}

		for(i=0;i<nbstates;i++) _out[i] /= sum;

	}

LFHTEMP	void ClassifierV<C>::EMinit(){
		int nbstates = classes.size();

		int i; for(i=0;i<nbstates;i++) ((ProbabilityDistribution<C>*)classes[i])->EMinit();
		if (unknown_scope != NULL){
			memset(unknown_scope+3,'\0',sizeof(double)*5);
		}

		if (prior_freq_learn != NULL){
			delete[](prior_freq_learn);
		}
	totdatabuffer =0;
		prior_freq_learn = new double[nbstates];
	if (fix_unknown_dens) unknown_scope[2] =0.0f;
	memset(prior_freq_learn, '\0', sizeof(double) *nbstates );
}
LFHTEMP	void ClassifierV<C>::EMAlphainit(double alp){
		int nbstates = classes.size();

		int i; for(i=0;i<nbstates;i++) ((ProbabilityDistribution<C>*)classes[i])->EMAlphainit(alp);
		if (unknown_scope != NULL){
			memset(unknown_scope+3,'\0',sizeof(double)*5);
		}


	if (prior_freq_learn != NULL){
		delete[](prior_freq_learn);
	}
	totdatabuffer =0;
	prior_freq_learn = new double[nbstates];
	if (fix_unknown_dens) unknown_scope[2] =0.0f;

	memset(prior_freq_learn, '\0', sizeof(double) *nbstates );
	}

	LFHTEMP	unsigned int ClassifierV<C>::NB_prestep() const{
		unsigned int fout = 0;
		unsigned int i,tmp;
		for(i=0;i< classes.size();i++) {tmp = classes[i]->NB_prestep(); if (tmp > fout) fout = tmp;}
		return fout;
	}


LFHTEMP	void ClassifierV<C>::EMregist(C &instance, const double * prob, const double weight){
	int i;		int nbstates = classes.size();

		double pix[8];
		if (unknown_scope == NULL){
			for(i=0;i<nbstates;i++) ((ProbabilityDistribution<C>*)classes[i])->EMregist(instance, prob[i] * weight);

			for(i=0;i<nbstates;i++){
				prior_freq_learn[i] += weight * prob[i];
			}
			totdatabuffer += weight;

		}else if (fix_unknown_dens == true){
			// passive guess

			// factor out probability of artifact!
			double tmp,sum;
			sum =0.0f;
			for(i=0;i<nbstates;i++){ tmp = classes[i]->LL(instance) * prior_freq[i];
				if (tmp > 0.0f) sum += prob[i] * tmp;
			}


			sum = exp(sum);
			tmp = sum / (sum + unknown_scope[1] * unknown_scope[0]); // foreground prob

			for(i=0;i<nbstates;i++){
				prior_freq_learn[i] += tmp * weight * prob[i];
			}
			unknown_scope[2] += (1 - tmp) * weight;
			totdatabuffer += weight;

			for(i=0;i<nbstates;i++) ((ProbabilityDistribution<C>*)classes[i])->EMregist(instance, tmp * prob[i] * weight);

		}else if (unknown_scope[2] == 0.0f){

			// passive guess
			for(i=0;i<nbstates;i++) ((ProbabilityDistribution<C>*)classes[i])->EMregist(instance, prob[i] * weight);
			unknown_scope[3] += weight;

			for(i=0;i<nbstates;i++){
				prior_freq_learn[i] += weight * prob[i];
			}
			totdatabuffer += weight;

		}else {
			// the probability of unknown is only defined by the anterior

			double tmp,sum,tmp2;
			sum =0.0f;
			for(i=0;i<nbstates;i++){ tmp = classes[i]->LL(instance);
				if (tmp > 0.0f) sum += prob[i] * tmp;
			}


			sum = exp(sum);
			tmp = sum / (sum + unknown_scope[0]);

			for(i=0;i<nbstates;i++){
				prior_freq_learn[i] += tmp * weight * prob[i];
			}
			totdatabuffer += weight;






			/*
			 unknown_scope[4] += 1.0f / (1.0f + exp(3.0f * unknown_scope[2]) * sum);
			 unknown_scope[5] += 1.0f / (1.0f + exp(unknown_scope[2]) * sum);
			 unknown_scope[6] += 1.0f / (1.0f + exp(-unknown_scope[2]) * sum);
			 unknown_scope[7] += 1.0f / (1.0f + exp(-3.0f * unknown_scope[2]) * sum);
			 */
			if (unknown_scope[2] > 0.125f){


				for(i=0;i<nbstates;i++) ((ProbabilityDistribution<C>*)classes[i])->EMregist(instance, prob[i] * tmp * weight);
				tmp = sum / (sum + unknown_scope[0] * exp(-unknown_scope[2]));
				// far approach:	ignores likelihood function
				tmp2 = sum / (unknown_scope[0] * exp(-unknown_scope[2]));
				sum = unknown_scope[0] / sum;
				if ((!ExCo<double>::isValid(sum))||(!ExCo<double>::isValid(tmp2))){
					unknown_scope[3] += 1.0f;
					unknown_scope[4] += 2.0f;
					unknown_scope[6] += 2.0f;
				}else {

					tmp = 2.0f + (sum + tmp2) / cosh(3.0f*unknown_scope[2]);
					tmp2 = 2.0f + (sum + tmp2) / cosh(unknown_scope[2]);

					if ((ExCo<double>::isValid(tmp))&&(ExCo<double>::isValid(tmp2))){

						if (ExCo<double>::isValid(2.0f * sinh(unknown_scope[2]) / tmp2)){

							pix[0] = (1.0f + (sum / cosh(3.0f*unknown_scope[2]))) / tmp;
							pix[1] = tanh(3.0f*unknown_scope[2]) / tmp;
							pix[2] = (1.0f + (sum/cosh(unknown_scope[2]))) / tmp2;
							pix[3] = tanh(unknown_scope[2]) / tmp2;

							if ((ExCo<double>::isValid(pix[0]))&&(ExCo<double>::isValid(pix[1]))&&(ExCo<double>::isValid(pix[2]))&&(ExCo<double>::isValid(pix[3]))){
								unknown_scope[4] += pix[0]* weight; // X1 + X4
								unknown_scope[5] += pix[1]* weight; // X1 - X4 /
								unknown_scope[6] += pix[2]* weight; // X2 + X3 / 2
								unknown_scope[7] += pix[3]* weight; // X2 - X3 /

								unknown_scope[3] += weight;
							}
						}
					}}
			}else {
				//close approach:    Uses likelihood function

				for(i=0;i<nbstates;i++) ((ProbabilityDistribution<C>*)classes[i])->EMregist(instance, prob[i] * tmp);
				sum /= unknown_scope[0];

				tmp = 1.0f + sum;

				pix[0] = 1.0f / tmp;
				double tmp2 = 2.0f + sum + (1.0f / sum);
				pix[1] = 1.0f / tmp2;
				pix[2] = 1.0f / (tmp2 * tmp);
				pix[3] = 1.0f / (tmp2 * tmp * tmp);

				if ((ExCo<double>::isValid(pix[0]))&&(ExCo<double>::isValid(pix[1]))&&(ExCo<double>::isValid(pix[2]))&&(ExCo<double>::isValid(pix[3]))){

					unknown_scope[3] += weight;
					unknown_scope[4] += pix[0] * weight;

					unknown_scope[5] += pix[1] * weight;
					unknown_scope[6] += pix[2] * weight;
					unknown_scope[7] += pix[3] * weight;
				}
			}
		}
	}
LFHTEMP	void ClassifierV<C>::EMclear(){
	int nbstates = classes.size();
	unsigned int i;
	if (EMsclasses != NULL) {
		for(i=0;i<classes.size();i++) { delete(classes[i]);	classes[i] = EMsclasses[i];	}
		delete[](EMsclasses);
		EMsclasses = NULL;
	}

	for(i=0;i<nbstates;i++) classes[i]->EMclear();
	}

LFHTEMP	double ClassifierV<C>::EMLogLikelihood() const{
		return(0.0f);
	}

LFHTEMP	double ClassifierV<C>::EMfinit(){
		int i;
	double coef[4];
	int nbstates = classes.size();
	for(i=0;i<nbstates;i++) {prior_freq[i] = prior_freq_learn[i] / totdatabuffer;}
		if (unknown_scope != NULL){
			if (fix_unknown_dens == true){
				coef[0] = 0.0f; // the leading contant is ignored
				for(i=0;i<nbstates;i++) coef[0] += ((ProbabilityDistribution<C>*)classes[i])->EMfinit();
				unknown_scope[1] = unknown_scope[2] / totdatabuffer;
			} else if (unknown_scope[2] > 0.125f){
				relerr = log(unknown_scope[6] - unknown_scope[7]) - log(unknown_scope[3]) - log(unknown_scope[1]) -log(2);
				printf("Far Step: F[%e] = %e !\n", unknown_scope[0] * exp(-2*unknown_scope[2]), (unknown_scope[4] - unknown_scope[5]) / (unknown_scope[3]));
				printf("Far Step: F[%e] = %e !\n", unknown_scope[0], (unknown_scope[6] - unknown_scope[7]) / (unknown_scope[3]));
				printf("Far Step: F[%e] = %e !\n", unknown_scope[0] * exp(2*unknown_scope[2]), (unknown_scope[6] + unknown_scope[7]) / (unknown_scope[3]));
				printf("Far Step: F[%e] = %e !\n", unknown_scope[0] * exp(4*unknown_scope[2]), (unknown_scope[4] + unknown_scope[5]) / (unknown_scope[3]));


				if (relerr < log(4)) { // dont update if too much data is in unknown class
					coef[1] = unknown_scope[3] * unknown_scope[1];
					coef[2] = unknown_scope[6] - unknown_scope[7];
					coef[0] = unknown_scope[3] * log(unknown_scope[0]) + log(coef[1]*coef[1] - 0.0625f*coef[2]*coef[2]) - 2.0f * log(coef[1]); // likelihood of unknown, and of frequency
					for(i=0;i<nbstates;i++) coef[0] += ((ProbabilityDistribution<C>*)classes[i])->EMfinit();
				} else {
					printf("Unstable Far Step: F[] = %e !\n", (unknown_scope[6] - unknown_scope[7]) / (unknown_scope[3]));
					ExCo<double>::toMinimum(coef[0]); // minimum value for likelihood!

				}
			}else if (unknown_scope[2] == 0.0f){
				coef[0] = 0.0f; // the leading contant is ignored
				for(i=0;i<nbstates;i++) coef[0] += ((ProbabilityDistribution<C>*)classes[i])->EMfinit();
			}else{
				// close approach

				//	printf("Close Step:\n");
				relerr = log(unknown_scope[4]) - log(unknown_scope[1] * unknown_scope[3]);

				if ((relerr < log(4))&&(ExCo<double>::isValid(relerr))) { // dont update if too much data is in unknown class
					coef[1] = unknown_scope[3] * unknown_scope[1];
					coef[0] = unknown_scope[3] * log(unknown_scope[0]) + log(coef[1]*coef[1] - 0.25f*unknown_scope[4]*unknown_scope[4]) - 2.0f * log(coef[1]); // likelihood of unknown, and of frequency
					for(i=0;i<nbstates;i++) coef[0] += ((ProbabilityDistribution<C>*)classes[i])->EMfinit();
				}else {
					printf("Unstable Near Step: F[] = %e !\n", unknown_scope[4] / unknown_scope[3]);
					ExCo<double>::toMinimum(coef[0]);
				}
			}

		}else {coef[0] =0; for(i=0;i<nbstates;i++) coef[0] += ((ProbabilityDistribution<C>*)classes[i])->EMfinit(); }

		return(coef[0]);
	}

LFHTEMP	void ClassifierV<C>::EMaccept(bool do_accept){
    unsigned int i;
    int nbstates = classes.size();
    if (!do_accept){
        if (EMsclasses == NULL) {LFH_ALIVE; exit(1);} // should always accept first!
        for(i=0;i<nbstates;i++) delete(classes[i]);
    }else{
    double coef[4];



    if (EMsclasses == NULL) EMsclasses = new ProbabilityDistribution<C>*[nbstates];
    else for(i=0;i<nbstates;i++) delete(EMsclasses[i]);
    for(i=0;i<nbstates;i++) EMsclasses[i] = classes[i];


    if (unknown_scope != NULL){
        if (fix_unknown_dens == true){

        }else{
        unknown_scope[9] = unknown_scope[0];

        if (unknown_scope[2] > 0.125f){
            // far approach	ignores likelihood function!
            //			printf("Far Step:\n");
            //		printf("%e\t%e\t%e\t%e\n", (unknown_scope[4]-unknown_scope[5])/(unknown_scope[3]), (unknown_scope[6]-unknown_scope[7])/(unknown_scope[3]),(unknown_scope[6]+unknown_scope[7])/(unknown_scope[3]),(unknown_scope[4]+unknown_scope[5])/(unknown_scope[3]));

            coef[0] =  (-unknown_scope[4] + 9.0f * unknown_scope[6]) - 8.0f * unknown_scope[3] * unknown_scope[1];
            coef[1] =  ((unknown_scope[5]/-3.0f) + 9.0f * unknown_scope[7]) ;
            coef[2] =  (unknown_scope[4] - unknown_scope[6]);
            coef[3] =  ((unknown_scope[5]/3.0f) - unknown_scope[7]);
            //		printf("%e\t%e\t%e\t%e\n", coef[0], coef[1], coef[2], coef[3]);
            double shift = CubicRealRoot(coef, false);
            //printf("Far  Step: F[%e] = %e !\n",unknown_scope[0], unknown_scope[6]/(unknown_scope[3]));
            //		printf("%e\t%e\n", shift, unknown_scope[2]);

            if (!ExCo<double>::isValid(shift)){

                unknown_scope[8] = unknown_scope[0] * exp((unknown_scope[6] > unknown_scope[3]* unknown_scope[1] ? -2.0f : 4.0f)* unknown_scope[2]);
                unknown_scope[2] *= 4.0f;

            }else if (fabs(shift) > 3.0f){ // outside!
                unknown_scope[8] = unknown_scope[0] * exp((shift < 0 ? -2.0f : 4.0f)* unknown_scope[2]);
                unknown_scope[2] *= 4.0f;

            }else{
                unknown_scope[8] = unknown_scope[0] * exp( (shift +1.0f) * unknown_scope[2]);
                unknown_scope[2] *= fabs(unknown_scope[1] -(unknown_scope[5]+ unknown_scope[6])/(unknown_scope[3]));
            }


            relerr = log(unknown_scope[6] - unknown_scope[7]) - log(unknown_scope[3]) - log(unknown_scope[1]) -log(2);




        }else if (unknown_scope[2] == 0.0f){
            unknown_scope[2] = 20.0f;
            relerr = 20.0f;
            unknown_scope[8] = exp(coef[0] / unknown_scope[3]);
        }else{
            // close approach
            double tmp;
            //	printf("Close Step:\n");
            relerr = log(unknown_scope[4]) - log(unknown_scope[1] * unknown_scope[3]);
            if ((ExCo<double>::isValid(unknown_scope[6]))&&(ExCo<double>::isValid(unknown_scope[7]))&&(unknown_scope[6] != 0.0f)&&(unknown_scope[7] != 0.0f)&&(fabs(1.0f / unknown_scope[6]) != 0.0f) &&(fabs(1.0f / unknown_scope[7]) != 0.0f) ){



                coef[0] = unknown_scope[4] - unknown_scope[1] * unknown_scope[3];
                coef[1] = unknown_scope[5];
                coef[2] = unknown_scope[6];
                coef[3] = unknown_scope[7];
                //		poly[2] = 2 * statbuf[3] + statbuf[2];
                //		poly[3] = 6 * statbuf[4] + 6 * statbuf[3] + statbuf[2];
                //		printf("%e\t%e\t%e\t%e\n",coef[0], coef[1], coef[2] , coef[3]  );

                tmp = -CubicRealRoot(coef,false);

                //		printf("%e poly eval\n", coef[0] -tmp *(coef[1] -tmp * (coef[2] -tmp * coef[3])));

                if ((!ExCo<double>::isValid(tmp))||(fabs(tmp) > 2.0f)){

                    tmp = (unknown_scope[4] - unknown_scope[1] * unknown_scope[3]) / unknown_scope[5];
                    //				printf("%e\t%e\n",(unknown_scope[4] - unknown_scope[1] * unknown_scope[3]) ,  unknown_scope[5]);
                    if (!ExCo<double>::isValid(tmp)) tmp = 2.0f;

                }

            }else if (unknown_scope[3] > 0.0f){

                tmp = (unknown_scope[4] - unknown_scope[1] * unknown_scope[3]) / unknown_scope[5];
                //		printf("wrongwrong F[%e] = %e ! %e\n",tmp ,unknown_scope[3],unknown_scope[4]);
                //		printf("%e\t%e\n",(unknown_scope[4] - unknown_scope[1] * unknown_scope[3]) ,  unknown_scope[5]);
                if (!ExCo<double>::isValid(tmp)) tmp = 2.0f;


            } else {tmp = 2.0f; unknown_scope[8] = unknown_scope[0] * exp(-3); relerr = 4.0f;}



            if (fabs(tmp) >= 2.0f){
                unknown_scope[8] =  unknown_scope[0] * (tmp < 0.0f) ? exp(2.0f) : exp(-2.0f);
                unknown_scope[2] *= 1.25f;
                if (unknown_scope[2] > 20.0f) unknown_scope[2] = 0.0f;
            }else{
                unknown_scope[2] *= fabs(tmp) /2.0f;
                unknown_scope[8] = unknown_scope[0] * exp(-tmp);
            }


    }
    }
    }}}

LFHTEMP	void ClassifierV<C>::EMnext(double alpha){
    if (unknown_scope) {unknown_scope[0] = unknown_scope[9] * (1.0f - alpha) + unknown_scope[8] * alpha; }
    unsigned int i;
    if (EMsclasses != NULL) {
        for(i=0;i<classes.size();i++) classes[i] = EMsclasses[i]->EMnext(alpha);
    }
}
LFHTEMP	void ClassifierV<C>::update_to_defaut_transition(){
    unsigned int coor[2];
    unsigned int nbstates = classes.size();
    coor[1] = coor[0] = nbstates;
    transitions.setSizes(coor);

    if (prior_freq) {
        for(coor[0]=0;coor[0]<nbstates;coor[0]++)
            for(coor[1]=0;coor[1]<nbstates;coor[1]++) transitions(coor) = prior_freq[coor[0]];
    }else{
        for(coor[0]=0;coor[0]<nbstates;coor[0]++)
            for(coor[1]=0;coor[1]<nbstates;coor[1]++) transitions(coor) = ((double)1.0f) / nbstates;
    }
}
LFHTEMP	void ClassifierV<C>::setSequential(bool on_off){
    isseq = on_off;
    if (isseq) update_to_defaut_transition();
}
LFHTEMP	void ClassifierV<C>::setPrior(double* npr){
    unsigned int nbstates = classes.size();
    if (prior_freq) delete[](prior_freq);
    prior_freq = new double[nbstates];
    for(nbstates--;nbstates != 0xFFFFFFFF;nbstates--)  prior_freq = npr[nbstates];
    if (isseq){
        nbstates = classes.size();
        unsigned int coor[2];
        for(coor[0]=0;coor[0]<nbstates;coor[0]++)
            for(coor[1]=0;coor[1]<nbstates;coor[1]++) transitions(coor) = prior_freq[coor[0]];
    }
}
LFHTEMP	void ClassifierV<C>::push_back(ProbabilityDistribution<C>*  newdist) {
	classes.push_back(newdist);
	if (isseq){
		unsigned int coor[2];
		coor[1] = coor[0] =  classes.size();
		transitions.setSizes(coor);
	}

		unsigned int nbstates =classes.size() ;
	double* n = new double[nbstates];
	if (nbstates > 1) memcpy(n,prior_freq,sizeof(double)*(nbstates-1));

	delete[](prior_freq); prior_freq = n;
	prior_freq[nbstates-1] = 1.0f;
}
LFHTEMP	void ClassifierV<C>::set_prior(double* prfreq, double unkwown){
    unsigned int nbstates =classes.size();
    memcpy(prior_freq,prfreq,sizeof(double)*nbstates);
    if (fix_unknown_dens) unknown_scope[2] = unkwown;
}
LFHTEMP	ProbabilityDistribution<C>* ClassifierV<C>::pop_back(){
	unsigned int nbstates =classes.size()-1;
	ProbabilityDistribution<C>* fout = classes[nbstates];
	classes.pop_back();
	if (nbstates){
		if (isseq){
			unsigned int coor[2];
			coor[1] = coor[0] =  classes.size();
			transitions.setSizes(coor);
		}

		double* n = new double[nbstates];
		memcpy(n,prior_freq,sizeof(double)*nbstates);
		delete[](prior_freq);
		prior_freq = n;
	}
return fout;}


	LFHTEMP	double ClassifierV<C>::performEM(Vector<C> &obs, int nbstep, bool rand_init){
		// EM on independent observations

		double *likeli = new double[classes.size()];
		int i , j, x,r;
		double LLhood;
		double LLhood_old;
		double alpha = 1.0f;
		double sum,saf;

		printf("EM start! with %i data points!\n", obs.size());

		for(i=0;i<nbstep;i++){

			this->EMinit();
			if ((i==0)&&(rand_init)){
				for(x = 0;x< obs.size();x++) {
					r = rand() % classes.size();
					for(j=0;j< classes.size(); j++) likeli[j] = (r == j) ? 1.0f : 0.0f;
					this->EMregist(obs[x] , likeli);
				}
			}else{
			for(x = 0;x< obs.size();x++) {
				for(j=0;j< classes.size(); j++) likeli[j] = classes[j]->LL(obs[x]);

				saf = likeli[0];
				for(j=1;j< classes.size(); j++) if (saf < likeli[j]) saf = likeli[j];
				sum =0;
				for(j=0;j< classes.size(); j++) {likeli[j] = prior_freq[j] * exp(likeli[j] - saf); sum += likeli[j];}
				for(j=0;j< classes.size(); j++) likeli[j] /= sum;


				this->EMregist(obs[x] , likeli);
			}
			}
			LLhood = this->EMfinit();

			if (((LLhood_old > LLhood)||(!ExCo<double>::isValid(LLhood)))&&(i != 0)) {
				printf("Failed: Log-likelihood = %e\n", LLhood);
				alpha *= 0.1f;
				EMaccept(false);
			}else{
				printf("Step %i: Log-likelihood = %e\n", i , LLhood);
				LLhood_old = LLhood;
				alpha = 1.0f;
				EMaccept();

			}
			EMnext(alpha);
		}
		delete[](likeli);

		EMclear();

		return (LLhood);
	}


	LFHTEMP	template<int dims> void ClassifierV<C>::EMregist(DataGrid< C , dims> &data, DataGrid< double , dims >* weigths){
		unsigned int nbstates = classes.size();
		DataGrid< double, dims+1> prob;

		DataGrid< double, dims> unk_prob;

		Tuple<unsigned int, dims+1> coors;
		memcpy(&(coors[1]),data.dims, sizeof(unsigned int) * dims);
		coors[0] = nbstates;
		prob.setSizes(coors);
		unk_prob.setSizes(data.dims);

		unsigned int i;
		double* tmppt = new double[nbstates];
		double* off;

		// no HMM, simple loop and register!
		typename DataGrid< C, dims>::KeyIterator ite = data.getKeyIterator();

		if (ite.first()) do{
			double sum = 0.0f;
			ite.write_to(&(coors[1]));


			for(i=0;i<nbstates;i++) {classes[i](tmppt[i],data(ite()));
				if (!ExCo<double>::isValid(tmppt[i])) tmppt[i] = 0.0f;
				sum += tmppt[i]; prob(coors) = tmppt[i];
			}

			if (sum == 0.0f){
				for(coors[0]=0;coors[0]<nbstates;coors[0]++) prob(coors) = 1.0f / nbstates;
				unk_prob(ite()) = 0.0f;
			} else{
				for(coors[0]=0;coors[0]<nbstates;coors[0]++) prob(coors) = tmppt[coors[0]] / sum;
				if (unknown_scope) {
					if (weigths != NULL) unk_prob(ite()) = (*weigths)(ite()) / (1.0f + (unknown_scope[0] / sum) );
					else unk_prob(ite()) = 1.0f / (1.0f + (unknown_scope[0] / sum) ) ;
				}else unk_prob(ite()) = 1.0f;
			}
			for(i=0;i<nbstates;i++) tmppt[i] = off[i] / sum;

		}while(ite.next());


		memcpy(&(coors[1]),data.dims, sizeof(unsigned int) * dims);
		coors[0] = nbstates;

		switch(dims){
			case 1:{
				DataGrid< double, dims+1> back_buf1(coors);
				directionnal_HMM_learn(1, back_buf1, prob);
				}break;
			case 2:{
				DataGrid< double, dims+1> back_buf1(coors);
				DataGrid< double, dims+1> back_buf2(coors);

				directionnal_HMM_learn(1, back_buf1, prob);
				directionnal_HMM_learn(2, back_buf2, prob);

				directionnal_HMM_learn(2, prob, back_buf1);
				directionnal_HMM_learn(1, back_buf1, back_buf2);

				prob *= back_buf1;


				}break;
			case 3:{
				DataGrid< double, dims+1> back_buf1(coors);
				DataGrid< double, dims+1> back_buf2(coors);
				DataGrid< double, dims+1> back_buf3(coors);
				DataGrid< double, dims+1> back_buf4(coors);

				directionnal_HMM_learn(1, back_buf1, prob);
				directionnal_HMM_learn(2, back_buf2, prob);

				directionnal_HMM_learn(2, back_buf3, back_buf1);
				directionnal_HMM_learn(1, back_buf4, back_buf2);

				back_buf4 *= back_buf3;

				directionnal_HMM_learn(3, back_buf3, prob);

				directionnal_HMM_learn(3, prob, back_buf1);
				directionnal_HMM_learn(1, back_buf1, back_buf3);

				back_buf1 *= prob;

				directionnal_HMM_learn(3, prob, back_buf2);
				directionnal_HMM_learn(2, back_buf2, back_buf3);

				back_buf2 *= prob;

				directionnal_HMM_learn(3, back_buf3, back_buf4);
				directionnal_HMM_learn(2, back_buf4, back_buf1);
				directionnal_HMM_learn(1, prob, back_buf2);

				prob *= back_buf4;
				prob *= back_buf3;


				}break;
			default:
				printf("%i dim-HMM not supported here!!! (is ignored)\n", dims);
				swap = false;




		}


	//	if (swap) {

	//	if (ite.first()) do{
	//		ite.write_to(&(coors[1]));
	//		for(coors[0]=0;coors[0]<nbstates;coors[0]++) classes[coors[0]]->EMregist(data(ite()), back_buf1(coors) * unk_prob(ite()));
	//	}while(ite.next());

	//	} else {

			if (ite.first()) do{
				ite.write_to(&(coors[1]));
				for(coors[0]=0;coors[0]<nbstates;coors[0]++) classes[coors[0]]->EMregist(data(ite()), prob(coors) * unk_prob(ite()));
			}while(ite.next());

	//	}



		delete[](tmppt);

	}



	LFHTEMP	template<int dims> void ClassifierV<C>::ExpectHidden(DataGrid< double, dims+1> &prob, const DataGrid< C, dims> &data, bool include_unknown) const{
		/*

		 bool todel = (bypassHMM == NULL)&&(isseq);
		 if (!isseq) bypassHMM = &prob;
		 if (todel) bypassHMM = new DataGrid< double, dims+1>();
		 unsigned int nbstates = classes.size();

		 DataGrid< double, dims> unknown_map;
		 Tuple<unsigned int, dims+1> coors;
		 memcpy(& (coors[1]), &data.dims, sizeof(unsigned int) *dims);
		 coors[0] = nbstates;
		 bypassHMM->setSizes(coors);
		 if (unknown_scope != NULL) unknown_map->setSizes(&(coors[1]));
		 prob.setSizes(coors);
		 double* tmp = new double[nbstates];
		 typename DataGrid< C, dims>::KeyIterator ite = data.getKeyIterator();





		 if (ite.first()) do{
		 memcpy(& (coors[1]), & (ite()), sizeof(unsigned int) *dims);
		 sum = 0;
		 for(coors[0] = 0 ; coors[0]< nbstates;coors[0]++) {
		 (*(classes[coors[0]]))(tmp[coors[0]], data(ite()));
		 sum += tmp[coors[0]];
		 }

		 if (!isseq){
		 for(coors[0] = 0 ; coors[0]< nbstates;coors[0]++) (*bypassHMM)(coors) = tmp[coors[0]] / sum;
		 } else if (sum < ldexp(2.0f,-100)){
		 for(coors[0] = 0 ; coors[0]< nbstates;coors[0]++) (*(classes[coors[0]])) = tmp[coors[0]] * ldexp(2.0f,200);
		 }else if (sum > ldexp(2.0f,100)){
		 for(coors[0] = 0 ; coors[0]< nbstates;coors[0]++) (*(classes[coors[0]])) = tmp[coors[0]] * ldexp(2.0f,-200);
		 }else{
		 for(coors[0] = 0 ; coors[0]< nbstates;coors[0]++) (*(classes[coors[0]])) = tmp[coors[0]];
		 }

		 if (unknown_scope != NULL){
		 unknown_map(ite()) = 1.0f / ((sum / unknown_map[0]) + 1.0f);
		 }


		 } while (ite.next());

		 // marginals computed

		 if (isseq) {
		 // running HMM!

		 switch(dims){
		 case 1:

		 break;
		 case 2:






		 break;

		 }




		 }

		 if ((bypassHMM != NULL)&&(todel == false)){ // normalize





		 if (ite.first()) do{

		 memcpy(& (coors[1]), & (ite()), sizeof(unsigned int) *dims);
		 sum = 0;
		 for(coors[0] = 0 ; coors[0]< nbstates;coors[0]++) {
		 (*(classes[coors[0]]))(tmp[coors[0]], data(ite()));
		 sum += tmp[coors[0]];
		 }
		 if (sum < ldexp(2.0f,-100)){
		 for(coors[0] = 0 ; coors[0]< nbstates;coors[0]++) (*(classes[coors[0]])) = tmp[coors[0]] * ldexp(2.0f,200);
		 }else if (sum > ldexp(2.0f,100)){
		 for(coors[0] = 0 ; coors[0]< nbstates;coors[0]++) (*(classes[coors[0]])) = tmp[coors[0]] * ldexp(2.0f,-200);
		 }else{
		 for(coors[0] = 0 ; coors[0]< nbstates;coors[0]++) (*(classes[coors[0]])) = tmp[coors[0]];
		 }

		 if (unknown_scope != NULL){
		 unknown_map(ite()) = 1.0f / ((sum / unknown_map[0]) + 1.0f);
		 }

		 for(coors[0] = 0 ; coors[0]< nbstates;coors[0]++) (*bypassHMM)(coors) = tmp[coors[0]] / sum;

		 } while (ite.next());

		 }
		 // normalize and regist





		 delete[](tmp);
		 if (todel) delete(bypassHMM);*/
	}


LFHTEMP	template<int dims> void ClassifierV<C>::directionnal_HMM_learn(int dir, const DataGrid<double, dims> &marginals, DataGrid<double, dims> &posterior, DataGrid<double, 2> &tr){

	unsigned int coors[dims];
	unsigned int precoors[dims];


	unsigned int i,j,nbd;
	nbd = marginals.dims[0];
	double* tmppt = new double[nbd+1];
	double tmpdouble;

	// backward pass
	for(i=0;i < dims;i++) {coors[i] = marginals.dims[i]-1; precoors[i] = (i == dir) ? coors[i] +1 : coors[i];}

	do {

		if (precoors[dir] != marginals.dims[dir])
		// {			for(coors[0]=0;coors[0] < nbd;coors[0]++) posterior(coors) = 1.0f; }else
		{
			// = transit * ( marginals(precoor) \cdot posterior(precoor))

			for(precoors[0]=0;precoors[0] < nbd;precoors[0]++) tmppt[precoors[0]] = marginals(precoors) * posterior(precoors);
			tmppt[nbd] =0.0f;
			for(precoors[0]=0;precoors[0] < nbd;precoors[0]++){
				posterior(precoors) = tmppt[0] * tr(precoors[0],0);
				for(coors[0]=1;coors[0] < nbd ;coors[0]++) posterior(precoors) += tmppt[coors[0]] * tr(precoors[0],coors[0]);
				tmppt[nbd] += posterior(precoors);
			}
			if (tmppt[nbd] < ldexp(2, -100)) {
				tr *= 16.0f;
				for(coors[0]=0;coors[0] < nbd;coors[0]++) posterior(coors) *= ldexp(2, 200);
			}else if (tmppt[nbd] > ldexp(2, 100)) {
				tr *= 0.0625f;
				for(coors[0]=0;coors[0] < nbd;coors[0]++) posterior(coors) *= ldexp(2, -200);
			}

		}


		for(i=1;i<dims;i++) if (coors[i] == 0) { coors[i] = marginals.dims[i]-1; precoors[i] = (i == dir) ? coors[i] +1 : coors[i]; } else {coors[i]--;precoors[i]--; break;}
	}while (i<dims);



	// forward pass
	// time-reversible matrix has:


	do {

		// get forward in tmp
		tmppt[nbd] = 0;
		if (coors[dir] == 0){
			for(precoors[0]=0;precoors[0] < nbd;precoors[0]++){
				coors[0]=0;
				tmppt[precoors[0]] = prior_freq[coors[0]] * tr(0,precoors[0]);
				for(coors[0]++;coors[0] < nbd;coors[0]++) tmppt[precoors[0]] += prior_freq[coors[0]] * tr(coors[0],precoors[0]);
				coors[0] = precoors[0];
				tmppt[coors[0]] *= marginals(coors);
				tmppt[nbd] += tmppt[precoors[0]];
			}

		} else{

			for(precoors[0]=0;precoors[0] < nbd;precoors[0]++){
				coors[0]=0;
				tmppt[precoors[0]] = posterior(coors) * tr(0,precoors[0]);
				for(coors[0]++;coors[0] < nbd;coors[0]++) tmppt[precoors[0]] += posterior(coors) * tr(coors[0],precoors[0]);
				coors[0] = precoors[0];
				tmppt[coors[0]] *= marginals(coors);
				tmppt[nbd] += tmppt[precoors[0]];
			}

		}

		// change scale


		// learn here


		// update!
		tmppt[nbd] = 0;
		if (precoors[dir] == marginals.dims[i]){
			for(coors[0]=0;coors[0] < nbd;coors[0]++) tmppt[nbd] += tmppt[precoors[0]];
		}else{
			for(coors[0]=0;coors[0] < nbd;coors[0]++) {tmpdouble = posterior(precoors); posterior(precoors) = tmppt[precoors[0]]; tmppt[precoors[0]] *= tmpdouble; tmppt[nbd] += tmppt[precoors[0]];}
		}

		for(coors[0]=0;coors[0] < nbd;coors[0]++) {
			precoors[0] = coors[0];
			posterior(coors) = tmppt[precoors[0]] / tmppt[nbd];
		}


		// [v1 v2] v2 = V1 * T /cdot back[coor[dir] ].data[z]

		// forward update F' <- (F*T) \cdot MARG
		// Prob update <- F' *

		for(precoors[0]=0;precoors[0] < nbd;precoors[0]++){
			coors[0]=0;
			tmppt[nbd] = posterior(coors) * tr(0,precoors[0]);
			for(coors[0]++;coors[0] < nbd;coors[0]++) tmppt[nbd] += posterior(coors) * tr(coors[0],precoors[0]);
			posterior(precoors) *= tmppt[nbd];
		}

			for(precoors[0]=0;precoors[0] < nbd;precoors[0]++) tmppt[precoors[0]] = marginals(precoors) * posterior(precoors);
			tmppt[nbd] =0.0f;
			for(coors[0]=0;coors[0] < nbd;coors[0]++){
				posterior(coors) = tmppt[0] * tr(coors[0],0);
				for(precoors[0]=1;precoors[0] < nbd ;precoors[0]++) posterior(coors) += tmppt[precoors[0]] * tr(coors[0],precoors[0]);
				tmppt[nbd] += posterior(coors);
			}
			if (tmppt[nbd] < ldexp(2, -100)) {
				tr *= 16.0f;
				for(coors[0]=0;coors[0] < nbd;coors[0]++) posterior(coors) *= ldexp(2, 200);
			}else if (tmppt[nbd] > ldexp(2, 100)) {
				tr *= 0.0625f;
				for(coors[0]=0;coors[0] < nbd;coors[0]++) posterior(coors) *= ldexp(2, -200);
			}




		for(i=1;i<dims;i++) if (coors[i] == marginals.dims[i]-1) { coors[i] = 0; precoors[i] = (i == dir) ? 1 : 0; } else {coors[i]++;precoors[i]++; break;}
	}while (i<dims);



	//  1) all eigen values are 1 (preserved by averaging)
	//  2)



}

LFHTEMP	template<int dims> void ClassifierV<C>::directionnal_HMM_eval(int direction, const DataGrid<double, dims> &marginals, DataGrid<double, dims> &posterior, DataGrid<double, 2> &tr) const{



	/*
	// dim0 is for the states

	Tuple<double, nbstate> tmp_out;

	_out.setSizes(like.dims);

	Matrix<double,nbstate,nbstate> learn;

	int coor[nbdim];
	int z,w;

	Matrix<double,1,nbstate>* back = new Matrix<double,1,nbstate>[ like.dims[dir] ];
	Matrix<double,nbstate,1> tmpb[2];


	int swp=0;
	double tmp;
	int altdir;
	memset(coor,'\0',sizeof(unsigned int)*nbdim);


	Matrix<double,nbstate,nbstate> tmptransition = transition;

	//	memset(updatetransit,'\0',sizeof(double)*channels*channels); // asume it<s done

	while(true){

		coor[dir] = like.dims[dir];
		for(z=0;z<nbstate;z++)	back[coor[dir]-1].data[z] =1.0f;

		for(coor[dir]--;coor[dir]>0;coor[dir]--) {
			back[coor[dir]-1] = tmptransition * (back[coor[dir]].scale_rows( like(coor) ));
			tmp = back[coor[dir]-1].data[0];
			for(w=1;w<nbstate;w++) tmp += back[coor[dir]-1].data[w];
			if (tmp <= 0.000001f) tmptransition *= 2.0f;
			else if (tmp >= 100000.0f) tmptransition *= 0.5f;


		}


		for(z=0;z<nbstate;z++)	tmpb[swp].data[z] = boundary[z];

		for( ;coor[dir]<like.dims[dir];coor[dir]++) {
			Tuple<double, nbstate> curob = like(coor);
			tmpb[swp^1]= (tmpb[swp] * tmptransition).scale_cols( curob );

			tmp = tmpb[swp^1].data[0];
			for(w=1;w<nbstate;w++) tmp += tmpb[swp^1].data[w];
			if (tmp <= 0.000001f) tmptransition *= 2.0f;
			else if (tmp >= 100000.0f) tmptransition *= 0.5f;

			if ((coor[dir] == 0)||(dolearn == false)){
				for(z=0;z<nbstate;z++) tmp_out[z] = tmpb[swp^1].data[z] * back[coor[dir] ].data[z];
			}else{
				tmp =0;
				for(z=0;z<nbstate;z++) {
					tmp_out[z] = tmpb[swp^1].data[z] * back[coor[dir] ].data[z];


					for(w=0;w<nbstate;w++) {
						learn.data[w + z* nbstate] = tmpb[swp].data[z] * transition.data[w + z * nbstate] * curob[w] * back[coor[dir] ].data[w];
						tmp += learn.data[w + z* nbstate];
					}
				}
				if ((tmp != 0)&&(!isnan(tmp))) {learn *= (1.0f / tmp);

				learned_transition += WeightElem< Matrix<double,nbstate,nbstate>, 1 >(learn);}
			}
			swp ^=1;
			_out( coor ) =  tmp_out.normalize();

		}

		for(altdir=1;altdir<nbdim;altdir++) if (coor[(dir + altdir) % nbdim] >= like.dims[(dir + altdir) % nbdim]-1) coor[(dir + altdir) % nbdim] =0; else {coor[(dir + altdir) % nbdim]++; break;}
		if (altdir == nbdim) break;
	}



	*/

}


/*
#undef LFHTEMP
#define LFHTEMP template<class INPUT>


LFHTEMP DistanceError<INPUT>::DistanceError(): EMscope(NULL), inside_factor(0.0f){}

LFHTEMP void DistanceError<INPUT>::operator()(double & f_out , pair<INPUT, double> &instance) const{
	double error = ExOp::norm( (instance.first) - center) - radius; // deviation!
	// error is the expected value for
//	double lamba = pow(-0.5, 0.5 * nihvar); // std
//		lamba = 1.0f / ((error < lamba) ? lamba : error);


//	f_out = inside_factor + (1 - inside_factor) * erf( nihvar * error ); // probability for the
//	f_out *= inside_factor + lamba *exp(lamba * error);
	error = (error + instance.second) ; // deviation!
	f_out += fact *exp( error*nihvar*error);

}
LFHTEMP	void DistanceError<INPUT>::EMinit(){ // uses local gradient!
	if (!EMscope) EMscope = new pair<INPUT, Tuple<double,4> >();
	ExOp::zero(EMscope->first);
	ExOp::zero(EMscope->second);
}
LFHTEMP	void DistanceError<INPUT>::EMAlphainit(double alpha){
	if (!EMscope) {
		EMscope = new pair<INPUT, Tuple<double,4> >();
		ExOp::zero(EMscope->first);
		ExOp::zero(EMscope->second);
	}else{
		EMscope->first *= alpha;
		EMscope->second *= alpha;
	}
}
LFHTEMP	void DistanceError<INPUT>::EMregist(const pair<INPUT,double> &instance, const double prob){
	double norm = ExOp::norm( (instance.first) - center);

	double errorval  = norm + instance.second  ; // deviation!


	EMscope->second[0] += prob;
	EMscope->second[1] += prob * errorval;
	EMscope->second[2] += prob * errorval * errorval;

	errorval = nihvar * (errorval - radius) ; // deviation!

	EMscope->first += ((instance.first) - center) * (prob * errorval);



}
LFHTEMP	void DistanceError<INPUT>::EMclear(){
	if (EMscope) {delete(EMscope); EMscope = NULL;}
}

LFHTEMP	double DistanceError<INPUT>::EMLogLikelihood() const{
	return(0.0f);
}


LFHTEMP	void DistanceError<INPUT>::EMfinit(){



	radius = EMscope->second[1] / EMscope->second[0];
//	double nradius = radius - (EMscope->second[1] / EMscope->second[0]);

	nihvar = 0.5f * EMscope->second[0] /  (EMscope->second[2] - radius * EMscope->second[1]);

	double nvar = -0.5f * log(M_PI * 2.0f * nihvar); // log of (expected Likelihood)...

	if (EMscope->second[3] == 0.0f){
		EMscope->second[3] = nvar;
		nvar = 1.0f;
	}else{
		double tmp =EMscope->second[3] - nvar ;
		nvar = 1.0f - 0.99f * exp( -tmp*tmp );
	}

	//nvar = 1.0f - 0.99f * exp(- nihvar * nradius * nradius); // estimate of overlap of distributions, defines alpha for update!

	double step = ExOp::norm( EMscope->first );

	step = sqrt(nihvar) / step; // maxstepsize

	center += EMscope->first * step * nvar;

	// store likelihood, for next updates
}

LFHTEMP	void DistanceError<INPUT>::show(FILE* o) const{
	fprintf(o, "DistanceError Distribution: rad=%f +- %f\tcenter :[", radius, pow(-0.5f, 0.5f * nihvar));
	ExOp::show(center,o,2);
	fprintf(o, "]\n");
}
*/

/*
    generic optimization has the form argmax \theta F(x, \theta)
    the number of parameters can be large, but can be partitioned as a list
    due to discontinuities, some set of parameter might need special function updates.
*/
/*
class OptimizationScope{
public:
    // optimization is of the form \sum F(x,theta...) + R(theta1) + R(thera2) + ...

    // and gradien step needs to comupute current value and propose an update (can abort if value is worse)

    void run(ThreadBase &tb, myHashmap<string, Anything> param);
};*/

#include "stats.hpp"






