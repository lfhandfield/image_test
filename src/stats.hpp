/*
 * primstats_tem.hpp
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




	template<class C, unsigned int nbstates>
	void Classifier<C,nbstates>::setUnknownProbability(double val){
		if ((val == 0.0f)^(unknown_scope == NULL)){
			if (val == 0) {delete[](unknown_scope); unknown_scope = NULL;}
			else{
				unknown_scope = new double[8];
				unknown_scope[0] = 0.0f;
				unknown_scope[1] = val;
				unknown_scope[2] = 0.0f; // step_size
			}
		}
	}

	template<class C, unsigned int nbstates>
	double Classifier<C,nbstates>::UnknownProbability(C &instance, const Tuple<double, nbstates>& prob){
		int i;
		double sum =0.0f;
		double tmp;
		for(i=0;i<nbstates;i++){(*classes[i])(tmp,instance);
			if (tmp != 0.0f) sum += prob[i] * log(tmp);
		}
		return(unknown_scope[0] / (exp(sum) + unknown_scope[0]));
	}


	template<class C, unsigned int nbstates>
	void Classifier<C,nbstates>::operator()(Tuple<double, nbstates>& _out, C& _in )const{
		Tuple<double, nbstates> tmp;
		int i;
		for(i=0;i<nbstates;i++){
			(*classes[i])(tmp[i],_in);
			/*	if (!ExCo<double>::isValid(tmp[i])) {

			 ((GaussianDistribution<2> *)classes[0])->show(stdout);
			 ((GaussianDistribution<2> *)classes[1])->show(stdout);

			 printf("%f\t%f\n", _in[0], _in[1]);
			 printf("%f\n", tmp[i]);
			 printf("AN EARLY ERROR!\n");exit(1);
			 }*/
		}
		// tmp[nbstates] = unknown_scope == NULL ? 0.0f : unknown_scope[0];
		_out = tmp.normalize();


		for(i=0;i<nbstates;i++){


		}

	}

	template<class C, unsigned int nbstates>
	void Classifier<C,nbstates>::EMinit(){
		unsigned int i; for(i=0;i<nbstates;i++) ((ProbabilityDistribution<C>*)classes[i])->EMinit();
		if (unknown_scope != NULL){
			memset(unknown_scope+3,'\0',sizeof(double)*5);
		}
	}
	template<class C, unsigned int nbstates>
	void Classifier<C,nbstates>::EMAlphainit(double alp){
		unsigned int i; for(i=0;i<nbstates;i++) ((ProbabilityDistribution<C>*)classes[i])->EMAlphainit(alp);
		if (unknown_scope != NULL){
			memset(unknown_scope+3,'\0',sizeof(double)*5);
		}
	}
	template<class C, unsigned int nbstates>
	void Classifier<C,nbstates>::EMregist(C &instance, const Tuple<double, nbstates> prob, bool should_learn){
		int i;
		double pix[8];
		if (unknown_scope == NULL){
			if (should_learn) for(i=0;i<nbstates;i++) ((ProbabilityDistribution<C>*)classes[i])->EMregist(instance, prob[i]);
		}else{
			// the probability of unknown is only defined by the anterior
			double tmp,sum,tmp2;
			sum =0.0f;
			for(i=0;i<nbstates;i++){(*classes[i])(tmp,instance);
				if (tmp > 0.0f) sum += prob[i] * log(tmp);
			}

			sum = exp(sum);

			tmp = sum / (sum + unknown_scope[0]);





			/*
			 unknown_scope[4] += 1.0f / (1.0f + exp(3.0f * unknown_scope[2]) * sum);
			 unknown_scope[5] += 1.0f / (1.0f + exp(unknown_scope[2]) * sum);
			 unknown_scope[6] += 1.0f / (1.0f + exp(-unknown_scope[2]) * sum);
			 unknown_scope[7] += 1.0f / (1.0f + exp(-3.0f * unknown_scope[2]) * sum);
			 */
			if (unknown_scope[2] > 0.125f){
				if (should_learn) for(i=0;i<nbstates;i++) ((ProbabilityDistribution<C>*)classes[i])->EMregist(instance, prob[i] * tmp);
				// far approach
				tmp2 = sum / unknown_scope[0];
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
								unknown_scope[4] += pix[0];
								unknown_scope[5] += pix[1];
								unknown_scope[6] += pix[2];
								unknown_scope[7] += pix[3];
								unknown_scope[3] += 1.0f;
							}
						}
					}}
			}else if (unknown_scope[2] == 0.0f){
				// passive guess
				if (should_learn) for(i=0;i<nbstates;i++) ((ProbabilityDistribution<C>*)classes[i])->EMregist(instance, prob[i]);
				if (ExCo<double>::isValid(log(sum))){
					unknown_scope[3] += 1.0f;
					unknown_scope[4] += log(sum);

				}
			}else {
				//close approach
				if (should_learn) for(i=0;i<nbstates;i++) ((ProbabilityDistribution<C>*)classes[i])->EMregist(instance, prob[i] * tmp);
				sum /= unknown_scope[0];

				tmp = 1.0f + sum;

				pix[0] = 1.0f / tmp;
				double tmp2 = 2.0f + sum + (1.0f / sum);
				pix[1] = 1.0f / tmp2;
				pix[2] = 1.0f / (tmp2 * tmp);
				pix[3] = 1.0f / (tmp2 * tmp * tmp);

				if ((ExCo<double>::isValid(pix[0]))&&(ExCo<double>::isValid(pix[1]))&&(ExCo<double>::isValid(pix[2]))&&(ExCo<double>::isValid(pix[3]))){

					unknown_scope[3] += 1.0f;
					unknown_scope[4] += pix[0];

					unknown_scope[5] += pix[1];
					unknown_scope[6] += pix[2];
					unknown_scope[7] += pix[3];
				}
			}
		}
	}


	template<class C, unsigned int nbstates>
	void Classifier<C,nbstates>::EMclear(){
		int i; for(i=0;i<nbstates;i++) ((ProbabilityDistribution<C>*)classes[i])->clear();
	}

	template<class C, unsigned int nbstates>
	double Classifier<C,nbstates>::EMfinit(){
		int i;
		double LL_out = 0.0f;
		if (unknown_scope != NULL){
			double coef[4];

			if (unknown_scope[2] > 0.125f){
				// far approach
				//			printf("Far Step:\n");
				//		printf("%e\t%e\t%e\t%e\n", (unknown_scope[4]-unknown_scope[5])/(unknown_scope[3]), (unknown_scope[6]-unknown_scope[7])/(unknown_scope[3]),(unknown_scope[6]+unknown_scope[7])/(unknown_scope[3]),(unknown_scope[4]+unknown_scope[5])/(unknown_scope[3]));

				coef[0] =  (-unknown_scope[4] + 9.0f * unknown_scope[6]) - 8.0f * unknown_scope[3] * unknown_scope[1];
				coef[1] =  ((unknown_scope[5]/-3.0f) + 9.0f * unknown_scope[7]) ;
				coef[2] =  (unknown_scope[4] - unknown_scope[6]);
				coef[3] =  ((unknown_scope[5]/3.0f) - unknown_scope[7]);
				//		printf("%e\t%e\t%e\t%e\n", coef[0], coef[1], coef[2], coef[3]);
				double shift = CubicRealRoot(coef, false);
				//printf("Far Step: F[%e] = %e!\n",unknown_scope[0], unknown_scope[6]/(unknown_scope[3]));
				//		printf("%e\t%e\n", shift, unknown_scope[2]);



				shift *=antioss;
				if (!ExCo<double>::isValid(shift)){
					unknown_scope[0] *= exp((unknown_scope[6] > unknown_scope[3]* unknown_scope[1] ? -3.0f : 3.0f)* unknown_scope[2]);

				}else if (fabs(shift) > 3.0f){ // outside!
					unknown_scope[0] *= exp((shift < 0 ? -3.0f : 3.0f)* unknown_scope[2]);
					unknown_scope[2] *= 4.0f;

				}else{
					unknown_scope[0] *= exp(shift * unknown_scope[2]);
					unknown_scope[2] *= fabs(unknown_scope[1] -(unknown_scope[5]+ unknown_scope[6])/(unknown_scope[3]));
				}
				//		printf("end: %e\t    %e\n", unknown_scope[0], unknown_scope[2]);
				relerr = log(unknown_scope[6]) - log(unknown_scope[3] * unknown_scope[1]);

				if (relerr < 1.0f){ // dont update if too much data is in unknown class
					for(i=0;i<nbstates;i++) LL_out += ((ProbabilityDistribution<C>*)classes[i])->EMfinit();
					LL_out /= (1.0f - unknown_scope[1] * exp(relerr));
				} else {unknown_scope[0] *= exp(-3);antioss *=0.9f; }
			}else if (unknown_scope[2] == 0.0f){
				unknown_scope[0] = exp(unknown_scope[4] / unknown_scope[3] );
				unknown_scope[2] = 20.0f;
				relerr = 20.0f;
				antioss =1.0f;
				for(i=0;i<nbstates;i++) LL_out +=((ProbabilityDistribution<C>*)classes[i])->EMfinit();
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


				} else {tmp = 2.0f; unknown_scope[0] *= exp(-3);antioss *=0.9f; relerr = 4.0f;}
				//printf("Near Step: F[%e] = %e !\n",unknown_scope[0], unknown_scope[4] / unknown_scope[3]);



				tmp *=antioss;
				if (fabs(tmp) >= 2.0f){
					unknown_scope[0] *=  (tmp < 0.0f) ? exp(2.0f) : exp(-2.0f);
					unknown_scope[2] *= 1.25f;
					if (unknown_scope[2] > 20.0f) unknown_scope[2] = 0.0f;
				}else{
					unknown_scope[2] *= fabs(tmp) /2.0f;
					unknown_scope[0] *= exp(-tmp);
				}
				//		printf("step = %e, new guess %e\n", tmp,unknown_scope[0]);

				if ((relerr < 1.0f)&&(ExCo<double>::isValid(relerr))){ // dont update if too much data is in unknown class
					for(i=0;i<nbstates;i++) LL_out += ((ProbabilityDistribution<C>*)classes[i])->EMfinit();
					LL_out /= (1.0f - unknown_scope[1] * exp(relerr));
				}else {unknown_scope[0] *= exp(-3);antioss *=0.9f; } // want to converge from below
			}

		}else for(i=0;i<nbstates;i++) LL_out += ((ProbabilityDistribution<C>*)classes[i])->EMfinit();
		return(LL_out);
	}

	template<class C, unsigned int nbstates>
	Classifier<C,nbstates> Classifier<C,nbstates>::EMfinit(double alpha){
		int i;
		Classifier<C,nbstates> f_out;
		if (unknown_scope != NULL){
			double coef[4];

			if (unknown_scope[2] > 0.125f){
				// far approach
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



				shift *=antioss;
				if (!ExCo<double>::isValid(shift)){
					unknown_scope[0] *= exp((unknown_scope[6] > unknown_scope[3]* unknown_scope[1] ? -3.0f : 3.0f)* unknown_scope[2]);

				}else if (fabs(shift) > 3.0f){ // outside!
					unknown_scope[0] *= exp((shift < 0 ? -3.0f : 3.0f)* unknown_scope[2]);
					unknown_scope[2] *= 4.0f;

				}else{
					unknown_scope[0] *= exp(shift * unknown_scope[2]);
					unknown_scope[2] *= fabs(unknown_scope[1] -(unknown_scope[5]+ unknown_scope[6])/(unknown_scope[3]));
				}
				//		printf("end: %e\t    %e\n", unknown_scope[0], unknown_scope[2]);
				relerr = log(unknown_scope[6]) - log(unknown_scope[3] * unknown_scope[1]);

				if (relerr < 1.0f){ // dont update if too much data is in unknown class
					for(i=0;i<nbstates;i++) ((ProbabilityDistribution<C>*)classes[i])->EMfinit();

				} else {unknown_scope[0] *= exp(-3);antioss *=0.9f; }
			}else if (unknown_scope[2] == 0.0f){
				unknown_scope[0] = exp(unknown_scope[4] / unknown_scope[3] );
				unknown_scope[2] = 20.0f;
				relerr = 20.0f;
				antioss =1.0f;
				for(i=0;i<nbstates;i++) ((ProbabilityDistribution<C>*)classes[i])->EMfinit();
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


				} else {tmp = 2.0f; unknown_scope[0] *= exp(-3);antioss *=0.9f; relerr = 4.0f;}
				//printf("Near Step: F[%e] = %e !\n",unknown_scope[0], unknown_scope[4] / unknown_scope[3]);



				tmp *=antioss;
				if (fabs(tmp) >= 2.0f){
					unknown_scope[0] *=  (tmp < 0.0f) ? exp(2.0f) : exp(-2.0f);
					unknown_scope[2] *= 1.25f;
					if (unknown_scope[2] > 20.0f) unknown_scope[2] = 0.0f;
				}else{
					unknown_scope[2] *= fabs(tmp) /2.0f;
					unknown_scope[0] *= exp(-tmp);
				}
				//		printf("step = %e, new guess %e\n", tmp,unknown_scope[0]);

				if ((relerr < 1.0f)&&(ExCo<double>::isValid(relerr))){ // dont update if too much data is in unknown class
					for(i=0;i<nbstates;i++) ((ProbabilityDistribution<C>*)classes[i])->EMfinit();

				}else {unknown_scope[0] *= exp(-3);antioss *=0.9f; } // want to converge from below
			}

		}else for(i=0;i<nbstates;i++) f_out.push_back(ExOp::new_class(classes[i]->EMfinit(alpha)));
		return f_out;
	}


template<class C, unsigned int nbstates> void  Classifier<C,nbstates>::performEM(Vector< C > &obs, int nbstep){
	// EM on independent observations
	Tuple<double, nbstates> likeli;
	int i , j, x;
	double LLhood;
	double LLhood_old;
	double alpha = 1.0f;

	Vector< ProbabilityDistribution<C>* > swaps_classes;

	for(i=0;i<nbstep;i++){

		this->EMinit();
		for(x = 0;x< obs.size();x++) {
			for(j=0;j< classes; j++) classes[j](likeli[j], obs[i]);
			likeli.normalize();
			this->EMregist(obs[i] , likeli);
		}

		LLhood = this->EMfinit();

		if (LLhood_old > LLhood) {i--;
			printf("Failed: Log-likelihood = %f\n", LLhood);
			alpha *= 0.1f;
			for(j=0;j< classes; j++) delete(classes[j]);
		}else{
			printf("Step %i: Log-likelihood = %f\n", LLhood);
			LLhood_old = LLhood;
			alpha = 1.0f;
			if (swaps_classes.size() == nbstates){
			for(j=0;j< classes; j++) {delete(swaps_classes[j]); swaps_classes[j] =classes [j];}
			}else for(j=0;j< classes; j++) swaps_classes.push_back(classes [j]);
		}
		for(j=0;j< classes; j++) classes[j] = swaps_classes[j]->EMnext(alpha);
	}
}

#undef LFHTEMP
#define LFHTEMP template<unsigned int DIMS>

LFHTEMP EuclidianEllipse<DIMS>::EuclidianEllipse():fixvars(0), updateGreedy(false) {	alpha_search.init(1.0f, DIMS*2);}
LFHTEMP EuclidianEllipse<DIMS>::EuclidianEllipse(const Tuple<double, 2> &in_center): center(in_center), fixvars(0), updateGreedy(false){ alpha_search.init(1.0f, DIMS*2); ExOp::toZero(eccent);}
// LFHTEMP EuclidianEllipse<DIMS>::EuclidianEllipse(const EuclidianEllipse<DIMS>&o) : center(o.center),eccent(o.eccent),updateGreedy(o.updateGreedy), fixvars(o.fixvars),width(o.width),nihvar(o.nihvar),fact(o.fact),EMscope(o.EMscope),alpha_search(o.alpha_search){}

LFHTEMP void EuclidianEllipse<DIMS>::operator()(double & f_out , pair< Tuple<double, DIMS> , double> &instance) const{
	double error = ExOp::norm( (instance.first) - center + eccent) + ExOp::norm( (instance.first) - center - eccent) + instance.second - width;


	// error is the expected value for
	//	double lamba = pow(-0.5, 0.5 * nihvar); // std
	//		lamba = 1.0f / ((error < lamba) ? lamba : error);


	//	f_out = inside_factor + (1 - inside_factor) * erf( nihvar * error ); // probability for the
	//	f_out *= inside_factor + lamba *exp(lamba * error);

	f_out = fact *exp( error*nihvar*error);

}

LFHTEMP double EuclidianEllipse<DIMS>::LL(const pair< Tuple<double, DIMS> , double> &instance) const{
	double error = ExOp::norm( (instance.first) - center + eccent) + ExOp::norm( (instance.first) - center - eccent) + instance.second - width;
	return (log(fact) + error*nihvar*error);
}

LFHTEMP	void EuclidianEllipse<DIMS>::EMinit(){ // uses local gradient!
	ExOp::toZero(EMscope);
}
LFHTEMP	void EuclidianEllipse<DIMS>::EMAlphainit(double alpha){
	//		if (!EMscope) {
	//			EMscope = new pair<INPUT, Tuple<double,4> >();
	//			ExOp::zero(EMscope->first);
	//			ExOp::zero(EMscope->second);
	//		}else{
	//			EMscope->first *= alpha;
	//			EMscope->second *= alpha;
	//		}
}

LFHTEMP	void EuclidianEllipse<DIMS>::EMpreregist(unsigned int step, const pair< Tuple<double, DIMS> ,double> &instance, const double prob){
	// find center of mass!
	EMscope[0] += prob; // sw
	int i;
	for(i=0;i<DIMS;i++){
		EMscope[3 +i ] -= prob * instance.first[i];
	}

}

LFHTEMP void EuclidianEllipse<DIMS>::EMprefinit(unsigned int step){
	int i;
	for(i=0;i<DIMS;i++){
		center[i] = EMscope[3 +i] / EMscope[0];
		EMscope[3 +i] =0.0f;
	}
	EMscope[0] = 0.0f;
}

LFHTEMP	void EuclidianEllipse<DIMS>::EMregist(const pair< Tuple<double, DIMS> ,double> &instance, const double prob){


	Tuple<double, DIMS> tmppt = (instance.first) - center + eccent;
	double norma = tmppt.norm();
	tmppt = (instance.first) - center - eccent;
	double normb = tmppt.norm();

	double errorval = norma + normb + instance.second  ; // deviation!

	EMscope[2] += prob * errorval * errorval; // swx2
	errorval *= prob;
	EMscope[1] += errorval; // swx
	EMscope[0] += prob; // sw
	// sufficien stat for radius update, conditionn

	// likelifunc += sum prob * errorval /  )

	int i;
	for(i=0;i<DIMS;i++){
		EMscope[3 +i ] -= errorval * (((instance.first[i] - center[i] + eccent[i]) / norma) + ((instance.first[i] - center[i] - eccent[i]) / normb));
		EMscope[3 +i + DIMS] += errorval * (((instance.first[i] - center[i] + eccent[i]) / norma) - ((instance.first[i] - center[i] - eccent[i]) / normb));
		EMscope[3 +i + DIMS*2] -= prob * (((instance.first[i] - center[i] + eccent[i]) / norma) + ((instance.first[i] - center[i] - eccent[i]) / normb));
		EMscope[3 +i + DIMS*3] += prob * (((instance.first[i] - center[i] + eccent[i]) / norma) - ((instance.first[i] - center[i] - eccent[i]) / normb));
	}
}

LFHTEMP double EuclidianEllipse<DIMS>::EMfinit(){

	if (EMscope[0] == 0.0f) return(0.0f); // did not get any points, dont update

	if (updateGreedy){
		return 0.0f;



	}else{
	// updates radius and nihval only!
	width = EMscope[1] / EMscope[0];
	if (!this->has_fix_stddev()){
	nihvar = -0.5f * EMscope[0] / (EMscope[2] - width *  EMscope[1]);
	fact = sqrt(nihvar / -M_PI);
	}
	double derivative[DIMS*2];

	// cached next calculation
	int i;
	for(i=0;i<DIMS;i++){
		EMscope[3 + i] =        derivative[i]        = nihvar * (EMscope[3 + i]       - width *  EMscope[3 + i +DIMS*2]) / EMscope[0] ;
		EMscope[3 + i +DIMS]  = derivative[i + DIMS] = nihvar * (EMscope[3 + i +DIMS] - width *  EMscope[3 + i +DIMS*3]) / EMscope[0] ;
	}

	alpha_search.registAscent( (log(fact) - 0.5f) ,derivative);
	EMscope[2] = alpha_search(); // stores alpha
	//printf("wnabe alpha! %f\n", EMscope[2]);

	switch (fixvars & 3){
		case 0: return EMscope[0] * (log(fact) - 0.5f);
		case 2: return EMscope[0] * log(fact) - 0.5f * nihvar * (EMscope[2] - EMscope[1] * width); // log likelihood!
		default: return EMscope[0] * log(fact) - 0.5f * nihvar * (EMscope[2] - 2*EMscope[1] * width + EMscope[0] * width * width); // log likelihood!
		}
	}
}

LFHTEMP	double EuclidianEllipse<DIMS>::EMloglikelihood() const{
	switch (fixvars & 3){
		case 0: return EMscope[0] * (log(fact) - 0.5f);
		case 2: return EMscope[0] * log(fact) - 0.5f * nihvar * (EMscope[2] - EMscope[1] * width); // log likelihood!
		default: return EMscope[0] * log(fact) - 0.5f * nihvar * (EMscope[2] - 2*EMscope[1] * width + EMscope[0] * width * width); // log likelihood!
	}
}

LFHTEMP	EuclidianEllipse<DIMS> EuclidianEllipse<DIMS>::EMnext_exec(double alpha) const{
	// updates the focal points! alpha updates!

	EuclidianEllipse<DIMS> f_out;
	f_out.width = width;
	f_out.nihvar = nihvar;
	f_out.fact = fact;
	f_out.alpha_search = alpha_search;
	int i;
	for(i=0;i<DIMS;i++){
		f_out.center[i] = center[i] + alpha * (EMscope[2] * EMscope[3 + i] + ((0.01f * width * rand()) / RAND_MAX));
		f_out.eccent[i] = eccent[i] + alpha * (EMscope[2] * EMscope[3 + i + DIMS] + ((0.01f * width * rand()) / RAND_MAX)) ;
	}

	return f_out;
}

LFHTEMP	void EuclidianEllipse<DIMS>::EMclear(){
	//		if (EMscope) {delete(EMscope); EMscope = NULL;}
}

LFHTEMP	void EuclidianEllipse<DIMS>::show(FILE* o) const{
	fprintf(o, "DistanceError Distribution: rad=%f +- %f\tcenter :[", width, pow(-2.0f * nihvar, -0.5f ));
	ExOp::show(center,o,2);
	fprintf(o, "]+-[");
	ExOp::show(eccent,o,2);
	fprintf(o, "]\n");
}


#undef LFHTEMP
#define LFHTEMP template<unsigned int D>

LFHTEMP void DirichletDistribution<D>::updateKonstant(){
    double sum =param[0];
    konstant = lngamma(param[0]);
    for(int i=1;i<param.getSize();i++) {sum +=param[i]; konstant += lngamma(param[i]);}
    konstant -= lngamma(sum);
}
LFHTEMP DirichletDistribution<D>& DirichletDistribution<D>::toOne(){
    param.toOne();
    konstant = -lngamma((double)param.getSize());
    printf("konstant: %i, %e\n", param.getSize(), konstant);
return *this;}

LFHTEMP void DirichletDistribution<D>::operator()(double& _out, Tuple<double, D>& _in)const{
    unsigned int i;
    _out = konstant;
    for(i=0;i<param.getSize();i++) _out += param[i] * log(_in[i]);
}

LFHTEMP double DirichletDistribution<D>::getLL(const DataGrid<double, 2u> &samples)const{
    double fout = konstant * samples.dims[1];
    Tuple<uint32_t, 2u> coor;
    double tmp;
    for(coor[0]=0;coor[0]< samples.dims[0];coor[0]++){
        tmp =0;
        for(coor[1]=0;coor[1]< samples.dims[1];coor[1]++){
            tmp += log(samples(coor));
            if (!ExOp::isValid(log(samples(coor)))) printf("Err, f(%e) = %e\n", samples(coor),log(samples(coor)) );
        }
        fout += param[coor[0]] * tmp;
    }

return fout;}
LFHTEMP void DirichletDistribution<D>::wrLLDerivative(const DataGrid<double, 2u> &samples, Tuple<double, D> &fout)const{
    fout.setSize(samples.dims[0]);
    Tuple<uint32_t, 2u> coor;
    double alphasum = param[0];
    for(coor[0]++;coor[0]< samples.dims[0];coor[0]++) alphasum += param[coor[0]];
    alphasum = d_lngamma_dx(alphasum);
    for(coor[0]=0;coor[0]< samples.dims[0];coor[0]++){
        fout[coor[0]] = (d_lngamma_dx(param[coor[0]]) - alphasum) * samples.dims[1];
        for(coor[1]=0;coor[1]< samples.dims[1];coor[1]++){
            fout[coor[0]] += log(samples(coor));
        }
    }
}

/*
LFHTEMP void DirichletDistribution<D>::wrRandSample(Tuple< Tuple<double, DIMS> > & _fout) const{
    Tuple<double> gammafact;
    uint32_t i,j;
    for(i=0;i<param.getSize();i++) gammafact[i] = lngamma(param[i]);
    double sum;
    for(i=0;_fout.getSize();i++){
        sum =0;
        for(j=0;param.getSize();i++){
            _fout[i][j] =

            sum += _fout[i][j];
        }
        _fout[i] /= sum;
    }
}*/

#undef LFHTEMP
#define LFHTEMP template<unsigned int D>

LFHTEMP void ParabolicMixingDistribution<D>::setSize(int _size){
    minimum.setSize(_size);
    shape.setSize(_size);
}
LFHTEMP ParabolicMixingDistribution<D>& ParabolicMixingDistribution<D>::toOne(){
    minimum.toOne();
    shape.toZero();
    LL_constant = 0.0f;
    for(uint32_t i=0;i<minimum.getSize();i++) LL_constant += log((double)i);
return *this;}

LFHTEMP void ParabolicMixingDistribution<D>::operator()(double& _out, Tuple<double, D>& _in )const{
    _out = exp(LL_constant);
}
LFHTEMP double ParabolicMixingDistribution<D>::getLL(const DataGrid<double, 2u> &samples)const{
    Tuple<uint32_t, 2u> coor;
    //for(coor[1]=0;coor[1]< samples.dims[1];coor[1]++){


    //}
return LL_constant * samples.dims[1];}

LFHTEMP void ParabolicMixingDistribution<D>::updateKonstant(){

}

LFHTEMP template<Tuple_flag TF> double ParabolicMixingDistribution<D>::expectStep(KeyElem<uint32_t,Tuple<double, 0u, TF> > &output, const KeyElem<uint32_t,Tuple<double, 0u, TF> > &_in, const Tuple<double, 0u, TF> &ll_deriv, const Trianglix<double, 0u> &ll_dblder, const Tuple<double> &ll_jumps) const{ // update value and output LL difference
// this step needs to be guarantied to increase the likelihood! And it bounces around the hyper-tetrahedron boundary.
    unsigned int i;
    // Step 1: get direction for Newton step
    output = ll_dblder.mkInvDivi(ll_deriv);

    // Step 2: align direction within domain boundaries
    double sum = 1.0f;
    double dotprod = 0.0f;
    double bk = DBL_MIN;
    double fw = DBL_MAX;
    double tmp;
    unsigned int nbvalid =_in.getSize();
    for(i=0;i<_in.getSize();i++){
        if (_in[i] == 0.0f){
            nbvalid--;
            if (output[i] < 0.0f) output[i] =0.0f; // do not go in negative...
        }else {
            sum -= _in[i];
            tmp = _in[i] / output[i];
            if (!ExOp::isValid(tmp)) {output[i] = 0.0f; continue;}
            if (output[i] > 0.0f) {if (tmp > bk) bk = tmp;}
            else {if (tmp < fw) fw = tmp;}
        }
        dotprod += output[i];
    }
    tmp = sum / dotprod;
    if (dotprod > 0.0f){
        if ((sum <= 2.2204460492503131e-16 * _in.getSize())&&(dotprod > 0.0f)) { // make the projection in orthogonal direction 0
            for(i=0;i<_in.getSize();i++){
                if (_in[i] != 0.0f) output[i] -= dotprod / nbvalid;
            }
        }else if (tmp < fw) fw = tmp;
    }else if (tmp > bk) bk = tmp;
    // Step 3: line search in allowed direction, if any.
    dotprod = 0.0f;
    for(i=0;i<_in.getSize();i++) dotprod += output[i] * ll_deriv[i];
    sum = ll_dblder.Xformed_inner_product(output);

    if (sum >= 0.0f){ // Extrema are the best, find both values
        tmp = ( (sum * bk + dotprod) *bk > (sum * fw + dotprod) * fw ) ? bk : fw;
    }else{ // go toward minimum if within domain
        tmp = dotprod / sum;
        tmp = (tmp < bk) ? bk : (tmp > fw) ? fw : tmp;
    }
    for(i=0;i<_in.getSize();i++) output[i] = _in[i] + output[i] * tmp;
return (sum * tmp + dotprod) *tmp;}

#undef LFHTEMP
#define LFHTEMP template<class C, unsigned int R, unsigned int S>

LFHTEMP void BiClassifier<C,R,S>::setDistrib(Tuple<unsigned int, 2u> dims, Oper2<double, DataGrid<C, 2u> >* _distr){
    distrib = _distr;
    par_RM.setSizes(dims);
    par_D.setSizes(dims);
}

LFHTEMP void BiClassifier<C,R,S>::setRandom(int nbrows, int nbcols, int nbrowshidden, int nbcolumnshidden){
    Tuple<uint32_t, 2u> coor;
    if ((nbrows != 0)||(nbcols != 0)) {
        rowhid.setNBcols(nbrows);
        colhid.setNBcols(nbcols);
        par_colscale.setSize(nbcols).toZero();
        par_rowscale.setSize(nbrows).toZero();
        coor[0] = nbrowshidden;
        coor[1] = nbcolumnshidden;
        par_RM.setSizes(coor).toRand();
        par_D.setSizes(coor);
    } else{ printf("no hidden size!\n");exit(1);}
    double sum;

    // cheating here!
    /*for(coor[1]=0; coor[1] < rowhid.getNBcols(); coor[1]++){
        sum =10;
        for(coor[0] = 0; coor[0] < par_RM.sizes[0]; coor[0]++){
            sum += (rowhid(coor) = randExponential());
        }
        coor[0] = coor[1] / 40; rowhid(coor) += 10.0;
        for(coor[0]=0; coor[0] < par_RM.sizes[0]; coor[0]++) rowhid(coor) /= sum;
    }*/
    // cheating here!
    /*for(coor[1]=0; coor[1] < colhid.getNBcols(); coor[1]++){
        sum =10;
        for(coor[0] = 0; coor[0] < par_RM.sizes[1]; coor[0]++){
            sum += (colhid(coor) = randExponential());
        }
        coor[0] = coor[1] / 40; colhid(coor) += 10.0;
        for(coor[0]=0; coor[0] < par_RM.sizes[1]; coor[0]++) colhid(coor) /= sum;
    }*/

    for(coor[1]=0; coor[1] < rowhid.getNBcols(); coor[1]++){
        sum =0;
        for(coor[0] = 0; coor[0] < par_RM.sizes[0]; coor[0]++){
            sum += (rowhid(coor) = randExponential());
        }
        for(coor[0]=0; coor[0] < par_RM.sizes[0]; coor[0]++) rowhid(coor) /= sum;
    }
    for(coor[1]=0; coor[1] < colhid.getNBcols(); coor[1]++){
        sum =0;
        for(coor[0] = 0; coor[0] < par_RM.sizes[1]; coor[0]++){
            sum += (colhid(coor) = randExponential());
        }
        for(coor[0]=0; coor[0] < par_RM.sizes[1]; coor[0]++) colhid(coor) /= sum;
    }

    // uniform prior!
    /*
    rowhid_prior.setSize(nbrowshidden);
    if (auto ite = rowhid_prior.getIterator()) do{*ite = 1.0;} while(ite++);

    colhid_prior.setSize(nbcolumnshidden);
    if (auto ite = colhid_prior.getIterator()) do{*ite = 1.0;} while(ite++);*/


    // default parabolic prior!

    rowhid_prior_constant = 0;//-0.999999 / nbrowshidden;
    rowhid_prior.setSize(nbrowshidden);
    if (auto ite = rowhid_prior.getIterator()) do{*ite = (ite()[0] == (ite()[1])) ? 1.0: 0.0;} while(ite++);

    colhid_prior_constant = 0;//-0.999999 / nbcolumnshidden;
    colhid_prior.setSize(nbcolumnshidden);
    if (auto ite = colhid_prior.getIterator()) do{*ite = (ite()[0] == (ite()[1])) ? 1.0: 0.0;} while(ite++);
}


LFHTEMP void BiClassifier<C,R,S>::setSemiRandom(int nbrows, int nbrowshidden, Vector<uint32_t> input_col){
    Tuple<uint32_t, 2u> coor;
    int nbcols = input_col.getSize();
    coor[1]=0;
    for(coor[0]=0;coor[0]< input_col.getSize();coor[0]++) if (input_col[coor[0]] >= coor[1]) coor[1] = 1 + input_col[coor[0]];
    int nbcolumnshidden = coor[1];
    if ((nbrows != 0)||(nbcols != 0)) {
        rowhid.setNBcols(nbrows);
        colhid.setNBcols(nbcols);
        par_colscale.setSize(nbcols).toZero();
        par_rowscale.setSize(nbrows).toZero();
        coor[0] = nbrowshidden;
        coor[1] = nbcolumnshidden;
        par_RM.setSizes(coor).toRand();
        par_D.setSizes(coor);
    } else{ printf("no hidden size!\n");exit(1);}
    double sum;

    for(coor[1]=0; coor[1] < rowhid.getNBcols(); coor[1]++){
        sum =0;
        for(coor[0] = 0; coor[0] < par_RM.sizes[0]; coor[0]++){
            sum += (rowhid(coor) = randExponential());
        }
        for(coor[0]=0; coor[0] < par_RM.sizes[0]; coor[0]++) rowhid(coor) /= sum;
    }
    for(coor[1]=0; coor[1] < colhid.getNBcols(); coor[1]++) colhid.data[coor[1]][input_col[coor[1]]] = 1.0;

    // default parabolic prior!
    rowhid_prior_constant = 0;//-0.999999 / nbrowshidden;
    rowhid_prior.setSize(nbrowshidden);
    if (auto ite = rowhid_prior.getIterator()) do{*ite = (ite()[0] == (ite()[1])) ? 1.0: 0.0;} while(ite++);


    colhid_prior_constant = 0;//-0.999999 / nbcolumnshidden;
    colhid_prior.setSize(nbcolumnshidden);
    if (auto ite = colhid_prior.getIterator()) do{*ite = (ite()[0] == (ite()[1])) ? 1.0: 0.0;} while(ite++);
}

LFHTEMP ERRCODE BiClassifier<C,R,S>::initHidden(int nbrowshidden ,int nbcolumnshidden, int nbrows, int nbcols, double row_sparse_size, double col_sparse_size, bool reuseRows, bool reuseCols){
    if ((nbrowshidden <= 1)||(nbcolumnshidden <= 1)) {printf("no hidden size!\n");return 1;}
    if ((nbrows <= 1)||(nbcols <= 1)) {printf("no size for data!\n");return 1;}
    if (col_sparse_size < 1.0001f){printf("sparsity size is expected to be higher than 1.0\n"); return 1;}
    if (row_sparse_size < 1.0001f){printf("sparsity size is expected to be higher than 1.0\n"); return 1;}

    int i,j,k;
    Tuple<double> pri;
    Trianglix<double> tr;
    if (reuseRows){
        if (rowhid.data.getSize() != nbrows) {printf("unexpected length for hidden rows\n"); return 1;}
        row_states_prior.setDefaultRatesAndPrior(nbrowshidden, row_sparse_size -1.0);
        auto lscp = row_states_prior.mkEMscope();
        lscp.init();
        for(i=0;i < rowhid.getNBcols();i++) {
            if (rowhid.data[i].getSize() != 0) lscp.EMregist(rowhid.getColumn(i));
        }
        row_states_prior.learn(lscp,1.0);
    }else{
         pri.setSize(nbrowshidden).toOne();
         tr.setSize(nbrowshidden);
         for(k=0,j=0;j<nbrowshidden;j++){
             for(i=0;i<=j;i++) tr.data[k++] = ((double)(i +1))/(j +1);
         }
         row_states_prior.setRatesAndPrior(pri,tr,row_sparse_size -1.0);
         rowhid.data =  std::move(row_states_prior.mkSamples(nbrows));
    }

    if (reuseCols){
        if (colhid.data.getSize() != nbcols) {printf("unexpected length for hidden rows\n"); return 1;}
        col_states_prior.setDefaultRatesAndPrior(nbcolumnshidden, col_sparse_size -1.0);
        auto lscp = col_states_prior.mkEMscope();
        lscp.init();
        for(i=0;i < colhid.getNBcols();i++) {
            if (colhid.data[i].getSize() != 0) lscp.EMregist(colhid.getColumn(i));
        }
        col_states_prior.learn(lscp,1.0);
    }else{
        pri.setSize(nbcolumnshidden).toOne();
        tr.setSize(nbcolumnshidden);
        for(k=0,j=0;j<nbcolumnshidden;j++){
            for(i=0;i<=j;i++) tr.data[k++] = ((double)(i +1))/(j +1);
        }
        col_states_prior.setRatesAndPrior(pri,tr,col_sparse_size -1.0);
        colhid.data =  std::move(col_states_prior.mkSamples(nbcols));
    }

    Tuple<uint32_t, 2u> coor;
    coor[0] = row_states_prior.getSize();
    coor[1] = col_states_prior.getSize();
    par_RM.setSizes(coor).toRand();
    par_D.setSizes(coor);
    if (par_colscale.getSize() != nbcols) par_colscale.setSize(nbcols).toZero();
    if (par_rowscale.getSize() != nbrows) par_rowscale.setSize(nbrows).toZero();

return 0;}

LFHTEMP void BiClassifier<C,R,S>::run2D_EM(const DataGrid<C, 2u>& data){ // Not Working
    Tuple<uint32_t, 2u> coor;
    rowhid.setNBcols(data.dims[0]);
    colhid.setNBcols(data.dims[1]);
    coor[0] = par_RM.sizes[0];
    coor[1] = par_RM.sizes[1];
    ((MatrixDistribution*)distrib)->setSizes(coor, true);
    // Randomly sample hidden variables using default prior

    for(int iff=0; iff<10;iff++){


        this->setRandom(0,0);

    double* inds;
   // ((MatrixDistribution*)distrib)->EMinit();
  //  for(coor[0] = 0; coor[0] < colhid.dims[0]; coor[0]++){
   //     inds = colhid.data + coor[0] *  colhid.dims[0];
    //    for(coor[1]=0; coor[1] < rowhid.dims[1]; coor[1]++){
  //          ((MatrixDistribution*)distrib)->EMregistOuter( rowhid.data + coor[0] * rowhid.dims[0], inds   , data(coor), 1.0f);
   //     }
    //}



    ((MatrixDistribution*)distrib)->EMHiddenMatrix(data,  rowhid, colhid);

    /*
      2dEM:
      iteratively update either {Z,} {H}, or {parameter and metaparameters}
      Hence, probability distributions needs to know what to do in each step, and optionally cache data for the last step...


    */
    }
}
LFHTEMP void BiClassifier<C,R,S>::run2D_EM(const SparseMatrix<C>& data){ // predecated
    Tuple<uint32_t, 2u> coor;
    rowhid.setNBcols(data.getNBrows());
    colhid.setNBcols(data.getNBcols());
    coor[0] = par_RM.sizes[0];
    coor[1] = par_RM.sizes[1];
    ((MatrixDistribution*)distrib)->setSizes(coor, true);

    // Randomly sample hidden variables using default prior

    this->setRandom(0,0);
    ((MatrixDistribution*)distrib)->EMHiddenMatrix(data, rowhid, colhid);



}


LFHTEMP void BiClassifier<C,R,S>::initScale(const SparseMatrix<C>& data, bool doRows, bool doCols){
    uint32_t i,j,k;
    myHashmap<uint32_t, double> toopure;
    if (doCols){
        par_colscale.toZero();
        if (auto ite = data.getIterator()) do{
            if ((*ite) > 1){
                par_colscale[ite.getCol()][0] += 1.0;
            }
            par_colscale[ite.getCol()][1] += (*ite) ;
        }while (ite++);

        j=0;k=0;
        for(i=0;i<colhid.getNBcols();i++){
            if (par_colscale[i][0] == 0){ // only has ones... r and m would diverge, hidden state should be "pure" then
                if (par_colscale[i][1] == 0){
                    colhid.data[i].toMemfree();
                    k++;
                }else{
                    toopure[colhid.data[i].deref_key(rand() % colhid.data[i].getSize())] = 1.0;
                    colhid.data[i].toMemmove(toopure);
                    j++;
                }
                par_colscale[i][0] = nan("");
            }else{
                par_colscale[i][0] = 0;
                par_colscale[i][1] = 0;
            }
        }
        printf("%i half-unuseable cols, and %i empty ones\n",j,k);
    }
    if (doRows){
        par_rowscale.toZero();
        if (auto ite = data.getIterator()) do{
            if ((*ite) > 1){
                par_rowscale[ite.getRow()][0] += 1.0;
            }
            par_rowscale[ite.getRow()][1] += (*ite) ;
        }while (ite++);
        j=0;k=0;
        for(i=0;i<rowhid.getNBcols();i++){
            if (par_rowscale[i][0] == 0.0){ // only has ones... r and m would dverge
                if (par_rowscale[i][1] == 0.0){
                    rowhid.data[i].toMemfree();
                    k++;
                }else{
                    toopure[rowhid.data[i].deref_key(rand() % rowhid.data[i].getSize())] = 1.0;
                    rowhid.data[i].toMemmove(toopure);
                    j++;
                }
                par_rowscale[i][0] = nan("");
            }else{
                par_rowscale[i][0] = 0;
                par_rowscale[i][1] = 0;
            }
        }
        printf("%i half-unuseable rows, and %i empty ones\n",j,k);
    }
}

LFHTEMP template<class P> void BiClassifier<C,R,S>::run2D_EM_v2(const SparseMatrix<C>& data, const Vector<uint32_t> &excl_list,  P& tb, uint32_t nbstep, unsigned int nbthread){

 //   DataGrid<double, 2u> proj_X;
 //   Trianglix<double> row_w_cumul; row_w_cumul.setSize(par_RM.sizes[0]); // sum_i ZZ^t / Z^tUZ
 //   Trianglix<double> col_w_cumul; col_w_cumul.setSize(colprc.getSize()); // sum_j HH^t / H^tVH
   /* Tuple<double> meanproject;
    Tuple<double> rowsums[2];
    Tuple<double> colsums[2];
    TMatrix<double> Msums[2];
    Tuple<double> Msumstmp;
    Tuple<double> rowdenum[2];
    Tuple<double> coldenum[2];

    uint32_t cached_index;
    for(cached_index =0; cached_index <2u;cached_index++){
        rowsums[cached_index].setSize(data.getNBrows());
        colsums[cached_index].setSize(data.getNBcols());
        rowdenum[cached_index].setSize(data.getNBrows());
        coldenum[cached_index].setSize(data.getNBcols());
        Msums[cached_index].setSizes(par_RM.sizes[0],par_RM.sizes[1]);
    }cached_index =0;
    Msumstmp.setSize(par_RM.sizes[0]);*/

    tb.toSize(nbthread-1);
    Tuple<Parallel_EMSCOPE> parascope;
    Shared_EMSCOPE scp(*this,data,excl_list);

    scp.ignore_drop = false;
    scp.hiddenalpha = 1.0;

    double that_ll_constant=0.0;
    if (auto ite = data.getIterator()) do{
        that_ll_constant -= lngamma(1.0 + *ite);
    }while(ite++);

    Tuple<uint32_t, 2u> coor, coor2, coorC, coorR;
    coor[0] = par_RM.sizes[0];
    coor[1] = par_RM.sizes[1];
    TMatrix< Tuple<double ,2u> > oldRM;
    double old_LL = -DBL_MAX;

    uint32_t i,j,k;
    parascope.setSize(nbthread);
    for(i=0;i<nbthread;i++) {
        parascope[i].scp = &scp;
        parascope[i].setSizes(coor,0);
        //parascope[i].setPointers((Tuple<double>*)rowdenum,(Tuple<double>*)coldenum, &cached_index);
        parascope[i].setRange((i * scp.trdata.getNBcols()) / nbthread, ((i+1) * scp.trdata.getNBcols()) / nbthread,
                              (i * data.getNBcols()) / nbthread, ((i+1) * data.getNBcols()) / nbthread);
        parascope[i].alpha = 0.001;
        parascope[i].use_stdout = (i == nbthread-1);
    }
//    meanproject.setSize(par_RM.sizes[0]);

    //proj_X.setSizes(coor);
    //proj_X.toZero();


    //bss.initRowIndexFromHidden(rowhid);
    //bss.initColIndexFromHidden(colhid);
    //bss.initCounts(data, rowh[1],colh[1]);

    ProgressBarPrint pbar(20);

    double tmp,tmpsum;

    //DataGrid<double> dll_dZ_stat;
    //DataGrid<double> dll_dH_stat;


    //coor[0] = colhidden.dims[0]; coor[1] = data.getNBrows(); dll_dZ_stat.setSizes(coor);

    double ll,ll_best;

    FunctionScope<double&,double&,double> funabsoft(ExCo<double>::toAbsoftInverse,pow(0.5, 256.0f));
    Trianglix<double,2u> single_denum;
    Tuple<double, 2u> single_num;

    Tuple<double,0u,TUPLE_FLAG_REMOTE_MEMORY> dainput;

    uint32_t rowite;
    int step;
    CurvatureSearchScope cs;
    //scp.istimefor_columns = !scp.istimefor_columns;
    uint32_t fuz =0;
    double expected =0.0;
    double expected2 =0.0;
    double dropLL;
    double last_dll;
    double rectsize;



    Tuple<double,0u> toproj; toproj.setSize(2 * par_RM.sizes[0] * par_RM.sizes[1]);
    Tuple<double,0u> ddproj; ddproj.setSize(2 * par_RM.sizes[0] * par_RM.sizes[1]);
    cs.init(2 * par_RM.sizes[0] * par_RM.sizes[1], 1.0);


    double logalpha = -1.0;

    dropLL = 0.0f;
    if (auto ite = par_D.getIterator()) do{
        (*ite) = ((double)scp.maxgrid_nonzero_count(ite())) / (scp.nb_maxrowh[ite()[0]] * scp.nb_maxcolh[ite()[1]]); // p = c_ij / (s_i*t_j)
    }while(ite++);

    while(true){
        scp.newhidden.setNBcols((scp.istimefor_columns)? colhid.getNBcols() : rowhid.getNBcols());
        scp.newscale.setSize(scp.newhidden.getNBcols());
 //       scp.new_nbmax_h.setSizes((scp.istimefor_columns)? scp.row_nbmax_hcol.getDims() : scp.col_nbmax_hrow.getDims()  );

        scp.tmptmp = false; //(fuz >= (nbstep/2)); // heavy/light hidden prior
        //scp.heuristictime = true;

        if ((fuz % 10) == 5) doRecenter(); // hehe

        for(i=nbthread-1;i>0;i--) tb.startEvent(&(parascope[i]));
        tb.startEvent_ThenWait(&(parascope[i]));
        for(i=1;i<nbthread;i++) {
            ExOp::toAdd(parascope[0].ll,parascope[i].ll);
            parascope[0].alpha += parascope[i].alpha;
            parascope[0].nonzero_diff += parascope[i].nonzero_diff;
            parascope[0].nbmax_diff += parascope[i].nbmax_diff;
            parascope[0].nb_fail += parascope[i].nb_fail;
            parascope[0].nb_insist += parascope[i].nb_insist;
        }
        if (!scp.ignore_drop) {
            parascope[0].nonzero_diff += scp.maxgrid_nonzero_count;
            dropLL =0.0;
            if (scp.istimefor_columns){
                for(i=0;i<par_RM.sizes[1];i++) parascope[0].nbmax_diff[i] += scp.nb_maxcolh[i];
                if (auto ite = par_D.getIterator()) do{
                    rectsize = ((scp.nb_maxrowh[ite()[0]] * parascope[0].nbmax_diff[ite()[1]]));
                    if ((parascope[0].nonzero_diff(ite()) != 0)&&(parascope[0].nonzero_diff(ite()) != scp.nb_maxrowh[ite()[0]] * parascope[0].nbmax_diff[ite()[1]])){
                        dropLL += log((double)parascope[0].nonzero_diff(ite())) * parascope[0].nonzero_diff(ite());
                        dropLL += log((double)(rectsize - parascope[0].nonzero_diff(ite()))) * (rectsize - parascope[0].nonzero_diff(ite()));
                        dropLL -= log((double)(rectsize)) * rectsize;
                        dropLL += lngamma((double)rectsize);
                        dropLL -= lngamma((double)(rectsize- parascope[0].nonzero_diff(ite())));
                        dropLL -= lngamma((double)(parascope[0].nonzero_diff(ite())));
                    }// else iszero
                }while(ite++);

            }else{
                for(i=0;i<par_RM.sizes[0];i++) parascope[0].nbmax_diff[i] += scp.nb_maxrowh[i];
                if (auto ite = par_D.getIterator()) do{
                    rectsize = ((double)(parascope[0].nbmax_diff[ite()[0]] * scp.nb_maxcolh[ite()[1]]));
                    if ((parascope[0].nonzero_diff(ite()) != 0)&&(parascope[0].nonzero_diff(ite()) != parascope[0].nbmax_diff[ite()[0]] * scp.nb_maxcolh[ite()[1]])){
                        dropLL += log((double)parascope[0].nonzero_diff(ite())) * parascope[0].nonzero_diff(ite());
                        dropLL += log((double)(rectsize - parascope[0].nonzero_diff(ite()))) * (rectsize - parascope[0].nonzero_diff(ite()));
                        dropLL -= log(rectsize) * rectsize;
                        dropLL += lngamma(rectsize);
                        dropLL -= lngamma((double)(rectsize- parascope[0].nonzero_diff(ite())));
                        dropLL -= lngamma((double)(parascope[0].nonzero_diff(ite())));
                    }// else iszero
                }while(ite++);
            }
            //printf("dropoutLL %e -> %e\n", last_dll, dropLL);
            last_dll = dropLL;
            parascope[0].ll[0] += dropLL; // not exact (dropout), but not important
            parascope[0].ll[1] += dropLL;
        }
        parascope[0].alpha  /= nbthread;
        for(i=1;i<nbthread;i++) parascope[i].alpha = parascope[0].alpha ;

        rectsize = ((double)parascope[0].nb_insist) / ((scp.istimefor_columns) ? scp.data.getNBcols() : scp.trdata.getNBcols());
        if (rectsize > 0.5) scp.hiddenalpha *= exp(-1);
        else scp.hiddenalpha = exp( 0.75 * log(scp.hiddenalpha) -0.25 * rectsize);
        rectsize = ((scp.istimefor_columns) ? scp.data.getNBcols() : scp.trdata.getNBcols());

        printf("%i Likelihood %e -> %e (%c In rate%e FF %e)\n", fuz, parascope[0].ll[0] + that_ll_constant, parascope[0].ll[1] + that_ll_constant, scp.istimefor_columns ? 'C' : 'R', ((double)parascope[0].nb_insist) / rectsize, ((double)parascope[0].nb_fail) / rectsize);
        printf("expected change %e, got %e  (ddmag %e)\n", expected + expected2, parascope[0].ll[0] - old_LL, expected2);

        if ((old_LL < parascope[0].ll[1])||(fuz == 0)||((scp.tmptmp)&&(fuz == (nbstep/2)))){
            expected = 0.0;
            expected2 = 0.0;
            oldRM = par_RM;
            old_LL = parascope[0].ll[1];
            if (scp.istimefor_columns){
                for(i=0;i< colhid.getNBcols();i++){
                    if ((colhid.data[i].getSize() == 0)||(scp.newhidden.data[i].getSize() == 0)) continue;
                    if (colhid.data[i].deref_key(0) != scp.newhidden.data[i].deref_key(0)){
                        for(j = 0; j < data.data[i].getSize();j++){
                            coor[1] = data.data[i].deref_key(j);
                            coor[0] = colhid.data[i].deref_key(0);
                            scp.row_nbmax_hcol(coor)--;
                            coor[0] = scp.newhidden.data[i].deref_key(0);
                            scp.row_nbmax_hcol(coor)++;
                        }
                    }
                }
                for(i=0;i<par_RM.sizes[1];i++) scp.nb_maxcolh[i] = parascope[0].nbmax_diff[i];
                colhid.toMemmove(scp.newhidden);
                par_colscale.toMemmove(scp.newscale);
             }else{
                for(i=0;i< rowhid.getNBcols();i++){
                    if ((rowhid.data[i].getSize() == 0)||(scp.newhidden.data[i].getSize() == 0)) continue;
                    if (rowhid.data[i].deref_key(0) != scp.newhidden.data[i].deref_key(0)){
                        for(j = 0; j < scp.trdata.data[i].getSize();j++){
                            coor[1] = scp.trdata.data[i].deref_key(j);
                            coor[0] = rowhid.data[i].deref_key(0);
                            scp.col_nbmax_hrow(coor)--;
                            coor[0] = scp.newhidden.data[i].deref_key(0);
                            scp.col_nbmax_hrow(coor)++;
                        }
                    }
                }
                for(i=0;i<par_RM.sizes[0];i++) scp.nb_maxrowh[i] = parascope[0].nbmax_diff[i];
                rowhid.toMemmove(scp.newhidden);
                par_rowscale.toMemmove(scp.newscale);
            }

            for(i=1;i<nbthread;i++) {
                parascope[0].dd_massive += parascope[i].dd_massive;
                parascope[0].dll_dMR2 += parascope[i].dll_dMR2;
            }
            // uptate D, close close form hurray!
            scp.maxgrid_nonzero_count = parascope[0].nonzero_diff;
            if (auto ite = par_D.getIterator()) do{
                (*ite) = ((double)scp.maxgrid_nonzero_count(ite())) / (scp.nb_maxrowh[ite()[0]] * scp.nb_maxcolh[ite()[1]]);
            }while(ite++);
            if (fuz >= nbstep) break;

            // fill symmetries... yea
            if (auto ite = parascope[0].dd_massive.getIterator()) do{
                if ((ite()[0] % par_RM.sizes[0]) > (ite()[1] % par_RM.sizes[0])){
                    (*ite) = parascope[0].dd_massive.cell( ((ite()[1] % par_RM.sizes[0]) + ((ite()[0] / par_RM.sizes[0]) * par_RM.sizes[0])) , ((ite()[0] % par_RM.sizes[0]) + ((ite()[1] / par_RM.sizes[0]) * par_RM.sizes[0])) );
                }
            }while(ite++);
            if (auto ite = parascope[0].dd_massive.getIterator()) do{
                if ((ite()[0] % (par_RM.sizes[0] * par_RM.sizes[1])) > (ite()[1] % (par_RM.sizes[0] * par_RM.sizes[1]))){
                    (*ite) = parascope[0].dd_massive.cell(ite()[1] % (par_RM.sizes[0] * par_RM.sizes[1]), (par_RM.sizes[0] * par_RM.sizes[1]) + ite()[0] % (par_RM.sizes[0] * par_RM.sizes[1]));
                }
            }while(ite++);

            /*
            i=0;
            if (auto ite = parascope[0].dll_dMR2.getIterator()) do{
                toproj[i++] = (*ite)[0];
                toproj[i++] = (*ite)[1];
            }while(ite++);
            cs.updateAscent(parascope[0].ll[1] - dropLL, &(par_RM.data[0][0]), &(toproj[0]));
            */

            i=0;
            if (auto ite = parascope[0].dll_dMR2.getIterator()) do{
                toproj[i + par_RM.sizes[0] * par_RM.sizes[1]] = (*ite)[1];
                toproj[i++] = (*ite)[0];
            }while(ite++);
            parascope[0].dd_massive.toEigenTransform(funabsoft);
            toproj = parascope[0].dd_massive * toproj;


            #ifdef Rcpp_hpp
          //  if (scp.istimefor_columns){
            i=0;
                if (auto ite = par_RM.getIterator()) do{
                    expected += exp(logalpha) * toproj[i] * parascope[0].dll_dMR2(ite())[0];
                    expected += exp(logalpha) * toproj[i + par_RM.sizes[0] * par_RM.sizes[1]] * parascope[0].dll_dMR2(ite())[1];
                    (*ite)[1] += exp(logalpha) * toproj[i + par_RM.sizes[0] * par_RM.sizes[1]];
                    (*ite)[0] += exp(logalpha) * toproj[i++];
                }while(ite++);
           // }

            #endif

            /*

            if (auto ite = par_RM.getIterator()) do{
                single_num[0] = parascope[0].dll_dMR2(ite())[0];
                single_num[1] = parascope[0].dll_dMR2(ite())[1];
                single_denum.data[0] = parascope[0].dll_dMR2(ite())[2];
                single_denum.data[1] = parascope[0].dll_dMR2(ite())[4];
                single_denum.data[2] = parascope[0].dll_dMR2(ite())[3];
                single_denum.toEigenTransform(funabsoft);
                single_num = single_denum * single_num;

                if (ExOp::isValid(single_num)){
                    single_num *= 0.1;
                    //single_num[1] = (single_num[1] + single_num[0]);
                    //single_num[0] = single_num[1];
                    if ((*ite)[0] + single_num[0] < -10.0f) single_num[0] = -10.0 - (*ite)[0]; // dont let R too far...

                    (*ite) += single_num;
                    expected += single_num[0] * parascope[0].dll_dMR2(ite())[0];
                    expected += single_num[1] * parascope[0].dll_dMR2(ite())[1];
                    expected2 += single_num[0] * single_num[0] * parascope[0].dll_dMR2(ite())[2];
                    expected2 += single_num[1] * single_num[1] * parascope[0].dll_dMR2(ite())[3];
                    expected2 += 2.0 * single_num[0] * single_num[1] * parascope[0].dll_dMR2(ite())[4];

                }else{
                   printf("produced nans from:\n");
                   parascope[0].dll_dMR2(ite()).show();
                }

            }while(ite++);*/
            //expected += expected2;
            scp.istimefor_columns = !scp.istimefor_columns;
            logalpha = (logalpha * 127.0) / 128.0;

        }else{ // failed... fall backward and do not update
            if (fuz >= nbstep + 10) break;
            printf("was tried:\n");
            par_RM.show();
            printf("backup: %e\n", logalpha);
            oldRM.show();
            printf("FAILFAILFAILFAILFAILFAILFAILFAILFAILFAILFAIL %e\n", old_LL + that_ll_constant); fflush(stdout);
            par_RM = (oldRM * 3.0 + par_RM) * 0.25;
            printf("FAILFAILFAILFAILFAILFAILFAILFAILFAILFAILFAIL %e\n", parascope[0].ll[0] + that_ll_constant); fflush(stdout);
            expected *= 0.25;
            logalpha -= 1.0;
        }
        fuz++;
     }



     if (scp.ignore_drop){ // update this, was not used
         scp.maxgrid_nonzero_count.toZero();
         if (auto ite = data.getIterator()) do{
            coor[0] = rowhid.data[ite()[0]].deref_key(0);
            coor[1] = colhid.data[ite()[1]].deref_key(0);
            scp.maxgrid_nonzero_count(coor)++;
         }while (ite++);
         scp.nb_maxrowh.toZero();
         scp.nb_maxcolh.toZero();
         for(coor[1]=0;coor[1]<rowhid.getNBcols();coor[1]++) {
            rowhid.makeMaximumAsFirst(coor[1]);
            scp.nb_maxrowh[rowhid.data[coor[1]].deref_key(0)]++;
        }
        for(coor[1]=0;coor[1]<colhid.getNBcols();coor[1]++){
            colhid.makeMaximumAsFirst(coor[1]);
            scp.nb_maxcolh[colhid.data[coor[1]].deref_key(0)]++;
        }
         if (auto ite = par_D.getIterator()) do{
            (*ite) = ((double)scp.maxgrid_nonzero_count(ite())) / (scp.nb_maxrowh[ite()[0]] * scp.nb_maxcolh[ite()[1]]);
         }while(ite++);
    }


    // Infer "density' of technical "d"ropouts

/*
    for(step=0;step<10;step++){
        // step 1: compute likelihood (only, but cache some calculation if pertinent for derivatives)
        ll = 0.0;
        dll_dZ_stat.toZero(); // dll_dH_stat.toZero();
        dainput.setSize(rowhidden.dims[0]);
        rowsums[cached_index].toZero();
        colsums[cached_index].toZero();
        Msums[cached_index].toZero();
        for(coor[1]=0;coor[1]<rowhidden.dims[1];coor[1]++){
            if (data.getRowSize(coor[1]) == 0) continue; // nothing to see...
            coor[0]=0; dainput = rowh[cached_index^1](coor); // set address to tuple
           // for(coor[0]=0;coor[0]<rowhidden.dims[0];coor[0]++) dainput[coor[0] ] = rowhidden(coor);
            rowdenum[cached_index][coor[1]] = rowprc.Xformed_inner_product(dainput);
            ll -= log(rowdenum[cached_index][coor[1]]) * data.getRowSize(coor[1]);
            //row_w_cumul += Trianglix<double>(dainput) / rowdenum[coor[1]];
        }
        dainput.setSize(colhidden.dims[0]);
        pbar.start("Compute E-step");
        for(coor[1]=0;coor[1]<data.getNBcols();coor[1]++){ pbar.update(coor[1],data.getNBcols());// column
            if (data.getColSize(coor[1]) == 0) continue; // nothing to see...
            coorC[1] = coor[1];
            coor[0]=0; dainput = colh[cached_index^1](coor); // set address to tuple
            //for(coor[0]=0;coor[0]<colhidden.dims[0];coor[0]++) dainput[coor[0] ] = rowhidden(coor);
            coldenum[cached_index][coor[1]] = colprc.Xformed_inner_product(dainput);
            ll -= log(coldenum[cached_index][coor[1]]) * data.getColSize(coor[1]);
            if (!ExOp::isValid(log(coldenum[cached_index][coor[1]]))) {printf("%e * %i\n",coldenum[cached_index][coor[1]], data.getColSize(coor[1])); dainput.show(); colprc.show();}

            //col_w_cumul += Trianglix<double>(dainput) / coldenum;
            for(coorR[0]=0;coorR[0]<rowprc.getSize();coorR[0]++){
                coorC[0]=0;
                meanproject[ coorR[0] ] = means.data[coorR[0] ] * colh[cached_index^1](coorC);
                for(coorC[0]++;coorC[0]<colprc.getSize();coorC[0]++) meanproject[ coorR[0]  ] +=  means.data[coorR[0]+ coorC[0] * rowprc.getSize()]  *colhidden(coorC);
            }
            Msumstmp.toZero();
            for(rowite=0;rowite<data.data[coor[1]].getSize();rowite++){ // row
                coor[0] = data.data[coor[1]].deref_key(rowite);
                coorR[1] = coor[0];
                tmp = data.data[coor[1]].deref(rowite);
                for(coorR[0]=0;coorR[0]<rowprc.getSize();coorR[0]++){
                    tmp -= meanproject[coorR[0]] * rowh[cached_index^1](coorR);
                }
                for(coorR[0]=0,coorC[0]=0;coorC[0] < rowhidden.dims[0];coorR[0]++,coorC[0]++) {
                    // dll_dH_stat(coorC) += rowhidden(coorR) * tmp / rowdenum[cached_index][coorR[1]]; // strange strange
                    Msumstmp[ coorR[0] ] += rowh[cached_index^1](coorR) * tmp / rowdenum[cached_index][coorR[1]]; // strange strange indeed
                }
            }
            for(coorC[0]=0,i=0;coorC[0]<colhidden.dims[0];coorC[0]++){
                for(coorR[0]=0;coorR[0]<rowhidden.dims[0];coorR[0]++){
                    Msums[cached_index].data[i++] += colhidden(coorC) * Msumstmp[ coorR[0] ] / coldenum[cached_index][coorC[1]];
                }
            }

            // Maximize H_j!

            for(rowite=0;rowite<data.data[coor[1]].getSize();rowite++){ // row
                coor[0] = data.data[coor[1]].deref_key(rowite);
                coorR[1] = coor[0];
                tmp = data.data[coor[1]].deref(rowite);
                for(coorR[0]=0;coorR[0]<rowprc.getSize();coorR[0]++){
                    tmp -= meanproject[coorR[0]] * rowh[cached_index](coorR);
                }

                for(coorR[0]=0,coorC[0]=0;coorR[0] < colhidden.dims[0];coorR[0]++,coorC[0]++) dll_dZ_stat(coorR) += colhidden(coorC) * tmp / coldenum[cached_index][coorC[1]];

                tmp = 0.5 * tmp * tmp;
                ll -= tmp / (rowdenum[cached_index][coor[0]] * coldenum[cached_index][coor[1]]);
                rowsums[cached_index][coor[0]] += tmp / coldenum[cached_index][coor[1]];
                colsums[cached_index][coor[1]] += tmp / rowdenum[cached_index][coor[0]];
            }
        } pbar.finish();

        printf("LL[E_H%i] = %e\n", step, ll);

        // Maximize Z_i!

        for(coor[1]=0;coor[1]<rowhidden.dims[1];coor[1]++){


        }
        printf("LL[E_Z%i] = %e\n", step, ll);

        pbar.start("Compute M-step");

        for(coor[1]=0;coor[1]<data.getNBcols();coor[1]++){ pbar.update(coor[1],data.getNBcols());// column



        }pbar.finish();

        printf("current!\n");
        printf("LL[M%i] = %e\n", step, ll);

        if (alpha.needsUndo(ll, step == 0)){
            // undo all, BUT Z and H
            means -= Msums[cached_index^1] * (0.75f * alpha[0]);
            rowprc.toAddProduct(dll_dU, -0.75f * alpha[1]);
            colprc.toAddProduct(dll_dV, -0.75f * alpha[2]);
            alpha.registUndo(0.25f);
        }else{

            // compute remaining statistics (what were not useful for computing LL, but needed in dLL)
            // update to Z and H uses gradient, *must* improve LL at current A,B,M,U,V (no backtrack allowed)
            dainput.setSize(rowhidden.dims[0]);
            for(coor[1]=0;coor[1]<data.getNBcols();coor[1]++){
                tmpsum =0.0f;
                row_w_cumul.toZero();
                for(rowite=0;rowite<data.data[coor[1]].getSize();rowite++){ // row
                    coor[0] = data.data[coor[1]].deref_key(rowite);
                    tmp =0;
                    for(coorR[0]=0,coorR[1]=coor[0];coorR[0]<proj_X.dims[0];coorR[0]++) tmp += ExOp::mkSquare(rowhidden(coorR));
                    row_w_cumul += Trianglix<double>(dainput) / rowdenum[cached_index][coor[1]];
                }
                Trianglix<double> tmptri = row_w_cumul.mkOuterMult(means);
                tmptri += rowprc * data.data[coor[1]].getSize();
            }

            dll_dU.toZero();
            for(coor[1]=0;coor[1]<rowhidden.dims[1];coor[1]++){
                tmpsum = ( (rowsums[cached_index][coor[1]] / rowdenum[cached_index][coor[1]]) - data.getRowSize(coor[1])) / rowdenum[cached_index][coor[1]];
                for(i=0,k=0;i<rowhidden.dims[0];i++){
                    coor[0] = i;
                    tmp = rowhidden(coor) * tmpsum;
                    for(coor[0]=0;coor[0]<=i;coor[0]++) dll_dU.data[k++] += rowhidden(coor) * tmp;
                }
            }
            dll_dV.toZero();
            for(coor[1]=0;coor[1]<colhidden.dims[1];coor[1]++){
                tmpsum = ( (colsums[cached_index][coor[1]] / coldenum[cached_index][coor[1]]) - data.getColSize(coor[1]) ) /coldenum[cached_index][coor[1]];
                for(i=0,k=0;i<colhidden.dims[0];i++){
                    coor[0] = i;
                    tmp = colhidden(coor) * tmpsum;
                    for(coor[0]=0;coor[0]<=i;coor[0]++) dll_dV.data[k++] += colhidden(coor) * tmp;
                }
            }

            //printf("P dLL_dU = "); ExOp::show(dll_dU);
            //printf("P dLL_dV = "); ExOp::show(dll_dV);
            //printf("P dLL_dM = "); ExOp::show(Msums[cached_index]);
            //tmp = this->evalLL_routine(data,rowhidden,colhidden,means,rowprc,colprc);
            //Trianglix<double> testest; testest.setSize(dll_dU.getSize());
            //printf("reeval LL = %e\n", tmp);
            for(coor[0]=0,i=0;coor[0]< dll_dU.getSize() ;coor[0]++){
                for(coor[1]=0; coor[1]<=coor[0]; coor[1]++,i++){
                    rowprc.data[i] += (coor[1] == coor[0]) ?  0.000001f : 0.0000005f;
                    testest.data[i] = this->evalLL_routine(data,rowhidden,colhidden,means,rowprc,colprc) -tmp;
                    rowprc.data[i] -= (coor[1] == coor[0]) ?  0.000001f : 0.0000005f;
                }
            }
            testest *= 1000000.0f;
            printf("R dLL_dU = "); ExOp::show(testest);
            testest.setSize(dll_dV.getSize());
            for(coor[0]=0,i=0;coor[0]< dll_dV.getSize() ;coor[0]++){
                for(coor[1]=0; coor[1]<=coor[0]; coor[1]++,i++){
                    colprc.data[i] += (coor[1] == coor[0]) ?  0.000001f : 0.0000005f;
                    testest.data[i] = this->evalLL_routine(data,rowhidden,colhidden,means,rowprc,colprc) -tmp;
                    colprc.data[i] -= (coor[1] == coor[0]) ?  0.000001f : 0.0000005f;
                }
            }
            testest *= 1000000.0f;
            printf("R dLL_dV = "); ExOp::show(testest);



            means +=  Msums[cached_index] * alpha[0];
            rowprc.toAddProduct(dll_dU, alpha[1]);
            colprc.toAddProduct(dll_dV, alpha[2]);
            cached_index = 1 - cached_index;
        }
        break;
    }*/
    // centering RM and coef (sum of rows/cols is 0, and average(row factor) = average(col)

    doRecenter();
}

LFHTEMP class BiClassifier<C,R,S>::Task_run2D_EM_v4{
    public:
    const BiClassifier<C,R,S> &target;
    const SparseMatrix<C>& data;
    SparseMatrix<C> trdata;
    Tuple<uint32_t> colrange,rowrange;

    std::mutex mtex_priscope;

    TMatrix<double> row_freq_hcol; // (par_RM.sizes[1] x data.size[0])
    TMatrix<double> col_freq_hrow; // (par_RM.sizes[0] x data.size[1])
    TMatrix<double> nonzero_freq_grid;
    Tuple<double> freq_rowh; //(par_RM.sizes[0])
    Tuple<double> freq_colh; //(par_RM.sizes[1])


    double dropout_likelihood_factor;


    //DataGrid<uint32_t> row_nbmax_hcol; // (par_RM.sizes[1] x data.size[0])
    //DataGrid<uint32_t> col_nbmax_hrow; // (par_RM.sizes[0] x data.size[1])
    Tuple< Tuple<double ,2u> > center_shift;

    // shared output
    SparseMatrix<double> newhidden;
    Tuple< Tuple<double ,2u> > newscale;
    Tuple< char > update_type;
    bool use_heavy_prior;


    /*
    Vector<uint32_t> indexes_changes;

    bool use_stdout;
    double alpha;

    double ll[2];
    Tuple<double>* rowdata;
    Tuple<double>* coldata;
    uint32_t* cachind;
    uint32_t nb_fail;
    uint32_t nb_insist;


    Tuple<int32_t> nbmax_diff;

    //shared
// (par_RM.sizes[0])

    DataGrid<uint32_t> row_nbmax_hcol; // (par_RM.sizes[1] x data.size[0])
    DataGrid<uint32_t> col_nbmax_hrow; // (par_RM.sizes[0] x data.size[1])





    SparseMatrix<double> introns;
    SparseMatrix<double> data;





    //bool heuristictime;


    */
    /*

    data.setNBcols(_data.getNBcols());
    introns.setNBcols(_data.getNBcols());
    myHashmap<uint32_t, void> which;
    for(uint32_t i=0;i<excl_list.getSize();i++) which.addEntry(excl_list[i]);

    if (do_incl_introns){
        for(uint32_t i = 0 ; i< _data.getNBcols();i++){
            if (auto ite = _data.data[i].getIterator()) do{
                if (which.find(ite()) != 0xFFFFFFFF) continue;
                data.data[i][ite()] = (*ite)[0] + (*ite)[1];
                introns.data[i][ite()] = (*ite)[1];
            }while(ite++);
        }
    }else{
        for(uint32_t i = 0 ; i< _data.getNBcols();i++){
            if (auto ite = _data.data[i].getIterator()) do{
                if (which.find(ite()) != 0xFFFFFFFF) continue;
                data.data[i][ite()] = (*ite)[0];
                introns.data[i][ite()] = (*ite)[1];
            }while(ite++);
        }
    }

    printf("fixing %i and %i nbrows so they match\n", trdata.data.getSize(), target.rowhid.getNBcols());
    while(trdata.data.getSize() < _target.rowhid.getNBcols()) trdata.data.push_back();
    init_routine(_target);

    */

    // thread local;
    Tuple<Tuple<double, 5> > ll;
    Tuple<double> alphas;
    Tuple< Tuple<uint32_t, 3u> > nb_flags; // nb_fail, nb_insist, nb_change;

    Tuple< TMatrix<double> > nonzero_freq_diff;
    //Tuple< TMatrix<int32_t> > nonzero_diff_dbg;

    TMatrix< Tuple<double ,2u> > oldRM;

    Tuple< Tuple<double> > freq_diff;


    Tuple< TMatrix< Tuple<double,2u> > > dll_dMR2;
    Tuple< Trianglix< double > > dd_massive;




    // flags:
    bool ignore_drop;
    bool tmptmp;

    uint32_t taskID;

    Task_run2D_EM_v4(const BiClassifier<C,R,S> &_target, const SparseMatrix<C>& _data, uint32_t _nbthread, bool _use_heavy_prior):target(_target), data(_data), use_heavy_prior(_use_heavy_prior){
        trdata = _data.mkTranspose();
        rowrange.setSize(_nbthread+1);
        colrange.setSize(_nbthread+1);
        alphas.setSize(_nbthread).toOne() *= 0.001;

        if (target.colhid.data.getSize() != data.getNBcols()){printf("nb col mismatch %i, %i", target.colhid.data.getSize(), data.getNBcols()); exit(1);}
        if (target.rowhid.data.getSize() != trdata.getNBcols()){
            if (target.rowhid.data.getSize() < trdata.getNBcols()) {printf("warning nb row mismatch %i, %i", target.rowhid.data.getSize(), trdata.getNBcols()); exit(1);}
            while (trdata.getNBcols() < target.rowhid.data.getSize())  trdata.data.push_back();
        }

        dll_dMR2.setSize(_nbthread);
        dd_massive.setSize(_nbthread);
        nb_flags.setSize(_nbthread);
        //nb_insist.setSize(_nbthread);
        //nb_change.setSize(_nbthread);
        ll.setSize(_nbthread);
        nonzero_freq_diff.setSize(_nbthread);
        //nonzero_diff_dbg.setSize(_nbthread);
        //nbmax_diff.setSize(_nbthread);
        freq_diff.setSize(_nbthread);
        Tuple<uint32_t, 2u> coor;
        coor[0] = _target.par_RM.sizes[0];
        coor[1] = _target.par_RM.sizes[1];
        int i;
        for(i=0;i<_nbthread;i++) {
            dll_dMR2[i].setSizes(coor[0],coor[1]);
            dd_massive[i].setSize(coor[1]*coor[0]*2);

            nonzero_freq_diff[i].setSizes(coor);
            //nonzero_diff_dbg[i].setSizes(coor);
            //nbmax_diff[i].setSize((coor[0] > coor[1]) ? coor[0] : coor[1]);
            freq_diff[i].setSize((coor[0] > coor[1]) ? coor[0] : coor[1]);
            rowrange[i] = (i * trdata.getNBcols()) / _nbthread;
            colrange[i] = (i * data.getNBcols()) / _nbthread;
        }
        rowrange[i] = trdata.getNBcols();
        colrange[i] = data.getNBcols();

        printf("allocated for %i\n", rowrange[i] > colrange[i] ? rowrange[i] : colrange[i]); fflush(stdout);
        update_type.setSize(rowrange[i] > colrange[i] ? rowrange[i] : colrange[i]);


        coor[0] = _target.par_RM.sizes[0];
        coor[1] = data.getNBcols();
      //  col_nbmax_hrow.toSizes(coor).toZero();
        col_freq_hrow.setSizes(coor).toZero();

        coor[0] = _target.par_RM.sizes[1];
        coor[1] = trdata.getNBcols();
      //  row_nbmax_hcol.toSizes(coor).toZero();
        row_freq_hcol.setSizes(coor).toZero();

        coor[0] = _target.par_RM.sizes[0];
        coor[1] = _target.par_RM.sizes[1];

        nonzero_freq_grid.setSizes(coor).toZero();

        freq_rowh.setSize(_target.par_RM.sizes[0]).toZero();
        freq_colh.setSize(_target.par_RM.sizes[1]).toZero();

        center_shift.setSize( (_target.par_RM.sizes[0]) > _target.par_RM.sizes[1] ? _target.par_RM.sizes[0] : _target.par_RM.sizes[1] ).toZero();

        taskID=0;
    }

    uint32_t operator()(uint32_t threadID){



        if ((taskID & 0xFFE) == 0){
            //nonzero_diff_dbg[threadID].toZero();
            ExOp::toZero(ll[threadID]);
            dd_massive[threadID].toZero();
            dll_dMR2[threadID].toZero();

            nonzero_freq_diff[threadID].toZero();
            //nbmax_diff[threadID].toZero();
            freq_diff[threadID].toZero();
            ll[threadID].toZero();
            nb_flags[threadID].toZero();
        }


        Vector<uint32_t> indexes_changes;
        Tuple<uint32_t,2u> coor, acoor;
        uint32_t i,j,dir,nbest;
        Tuple<Tuple<double,2u> > curproj[4]; // R H_j C_j
        double tmp;
        Tuple< Tuple<double, 2u>, 0u> future_statistics[3];
        Trianglix< Tuple<double, 3u>, 0u> future_statistic_trix[3];
        Tuple<double, 3u> trix_input;
        SparseTuple< Tuple<double, 2u> > future_factor;

        SparseTuple<double> candidates[2];
        SparseTuple<double>  can_backup;
        Tuple<double, 3u> lld;
        Tuple<double, 4u> lls;
        Tuple<double, 2u> coeff, coeff_test;
        double sum[4];

        FunctionScope<double&,double&,double> funabsoft(ExCo<double>::toAbsoftInverse,pow(0.5, 256.0f));

        Tuple<double, 2u> d_scale[2];
        Trianglix<double, 2u> dd_scale;
        Tuple<Tuple<double,2u>, 0u> rever;

        double prior_llterm[3];
        double Rval,Mval,logexpmin,expRterm,prediction;
        union LocalDef{
            double evilterms[6];
            KeyElem<double, Tuple<double, 5u> > evilterm_alt;
            LocalDef(){}
        };
        LocalDef ld;

        Tuple<double> RMS_deriv;
        Tuple<Tuple<Tuple<double,2u> > > rever_buf;
        Tuple<Tuple<double, 6u> > dd_buf;

        bool haschange;
        int insist=0;
        int change=0;
        double accuracy =0.0;
        uint32_t nb_funky=0u;
        Tuple<double, 4u> tosolve_tmp[3];

        Tuple<double, 3u> tosolve_num2, solution;
        Trianglix<double, 3u> tosolve_den2, tosolve_den2_tmp;
        Tuple<double> drop_buf;

        double expexp, expexp2;
        bool checksafe;
        Vector<int> torem;
        uint32_t wtorem;
        Tuple<uint32_t, 2u> mergesuggestion;
        uint32_t toremdirection;
        bool maindoutcmp;
        switch(taskID){
            case 0:

            RMS_deriv.setSize(target.par_RM.sizes[1]);

            rever_buf.setSize(target.rowhid.getNBcols());
            dd_buf.setSize(target.rowhid.getNBcols());

            for(i=0;i<3;i++) {
                future_statistics[i].setSize(target.par_RM.sizes[0]);
                future_statistic_trix[i].setSize(target.par_RM.sizes[0]);
            }
            ExOp::toZero(sum);
            for(coor[1]=colrange[threadID];coor[1] < colrange[threadID+1];coor[1]++){ thrbase.updateProgress(threadID);
                //printf("a%ib%ic\n", threadID, coor[1]);
                const SparseTuple<double>& curcol = target.colhid.getColumn(coor[1]);
                if (auto ite =curcol.getIterator()) do {if (ite() >= target.par_RM.sizes[1]) {printf("IMPOSSIBLBLBLB!!!\n"); curcol.show(); exit(1);} }while(ite++);
                if (!ExOp::isValid(target.par_colscale[coor[1]][0])){
                    newscale[coor[1]] = target.par_colscale[coor[1]];

                    if (data.getColSize(coor[1]) == 0) {update_type[coor[1]] = 'S'; candidates[0] = target.colhid.getColumn(coor[1]); newhidden.memmoveColumn(coor[1], candidates[0]); continue;} // nothing to see...
                    // column only has 1s, cannot be used to update R or M, drop out only, pure state guarantied

                    // trying a random alternatives

                    if (curcol.getSize() != 1) {printf("half-unsuable col%i should have 1 dim only... got %i\n", coor[1], curcol.getSize()); exit(1);}
                    candidates[0].toMemfree();
                     candidates[1].toMemfree();
                     i = rand() % (target.par_RM.sizes[1] - 1);
                     candidates[0][i] = 1.0;
                     if (i >= curcol.deref_key(0)) i++;

                     lls[0] = target.col_states_prior(candidates[0]) * ((use_heavy_prior) ? data.data[coor[1]].getSize() : 1);
                     lls[2] = target.col_states_prior(curcol) * ((use_heavy_prior) ? data.data[coor[1]].getSize() : 1);
                    if (target.par_RM.sizes[1] > 2){
                         j = rand() % (target.par_RM.sizes[1] - 2);
                         if (j >= curcol.deref_key(0)) j++;
                         if (j >= i) j++;
                        candidates[1][j] = 1.0;
                        lls[1] = target.col_states_prior(candidates[1]) * ((use_heavy_prior) ? data.data[coor[1]].getSize() : 1);
                     }else j=i;

                     ExOp::toZero(lld);

                     acoor[1] = curcol.deref_key(0);
                    for(acoor[0] =0;acoor[0]< target.par_RM.sizes[0]; acoor[0]++){
                        coor[0] = acoor[0];
                        if (freq_rowh[coor[0]] <= col_freq_hrow(coor)- 0.001) thrbase.msgs.insert(to_string(coor[0]) + string("got negative!") + to_string((double)col_freq_hrow(coor)));
                        lld[2] += col_freq_hrow(coor) * log(target.par_D(acoor)) + (freq_rowh[coor[0]] - col_freq_hrow(coor)) * log(1.0 - target.par_D(acoor));
                    }
                     acoor[1] = candidates[0].deref_key(0);
                    for(acoor[0] =0;acoor[0]< target.par_RM.sizes[0]; acoor[0]++){
                        coor[0] = acoor[0];
                        if (freq_rowh[coor[0]] <= col_freq_hrow(coor)- 0.001) thrbase.msgs.insert(to_string(coor[0]) + string("got negative!") + to_string((double)col_freq_hrow(coor)));
                        lld[0] += col_freq_hrow(coor) * log(target.par_D(acoor)) + (freq_rowh[coor[0]] - col_freq_hrow(coor)) * log(1.0 - target.par_D(acoor));
                    }
                     acoor[1] = j;
                    for(acoor[0] =0;acoor[0]< target.par_RM.sizes[0]; acoor[0]++){
                        coor[0] = acoor[0];
                        if (freq_rowh[coor[0]] <= col_freq_hrow(coor)- 0.001) thrbase.msgs.insert(to_string(coor[0]) + string("got negative!") + to_string((double)col_freq_hrow(coor)));
                        lld[1] += col_freq_hrow(coor) * log(target.par_D(acoor)) + (freq_rowh[coor[0]] - col_freq_hrow(coor)) * log(1.0 - target.par_D(acoor));
                    }
                     lld[0] *= dropout_likelihood_factor;
                     lld[1] *= dropout_likelihood_factor;
                     lld[2] *= dropout_likelihood_factor;
                     j = ((candidates[1].getSize() == 0)||(lls[0] + lld[0] > lls[1] + lld[1])) ? 0 : 1;
                     j = (lls[j] + lld[j] > lls[2] + lld[2]) ? j : 2;
                     if (j == 2){
                        candidates[0] = curcol; newhidden.memmoveColumn(coor[1], candidates[0]);
                        update_type[coor[1]] = 'S';
                     }else{
                        newhidden.memmoveColumn(coor[1], candidates[j]);
                        update_type[coor[1]] = 'D';

                        freq_diff[threadID][curcol.deref_key(0)] -= 1.0;
                        freq_diff[threadID][candidates[j].deref_key(0)] += 1.0;

                        for(i = 0; i < data.data[coor[1]].getSize();i++){
                            tosolve_num2[0] = 1.0 / (target.rowhid.data[ data.data[coor[1]].deref_key(i)].getSize());
                            if (auto ite2 = target.rowhid.data[ data.data[coor[1]].deref_key(i)].getIterator()) do{
                                acoor[0] = ite2();
                                acoor[1] = curcol.deref_key(0);        nonzero_freq_diff[threadID](acoor) -= tosolve_num2[0];
                                acoor[1] = candidates[j].deref_key(0); nonzero_freq_diff[threadID](acoor) += tosolve_num2[0];
                            }while(ite2++);
                        }
                        nb_flags[threadID][2]++;
                     }
                    //ll[threadID][0] += 0;
                    //ll[threadID][1] += 0;
                    ll[threadID][2] += lld[2];
                    ll[threadID][3] += lld[j];
                    ll[threadID][0] += lls[2];
                    ll[threadID][1] += lls[j];
                    ll[threadID][4] += lls[2];
                    candidates[1].toMemfree();
                    candidates[0].toMemfree();

                    if (auto ite =newhidden.data[coor[1]].getIterator()) do {if (ite() >= target.par_RM.sizes[1]) {printf("IMPOSSdddIBLBLBLB!!!\n"); curcol.show(); exit(1);} }while(ite++);


                    // put in stone... for now
                    /*
                    RM_deriv.toZero();

                    for(coor[0]=0;coor[0]<col_nbmax_hrow.dims[0];coor[0]++){
                        acoor[0] = coor[0];
                        for(acoor[1] =0;acoor[1]< target.par_RM.sizes[1]; acoor[1]++){
                            RM_deriv[acoor[1]][0] += col_nbmax_hrow(coor) * log(target.par_D(acoor));
                            RM_deriv[acoor[1]][0] += (nb_maxrowh[coor[0]] - col_nbmax_hrow(coor)) * log(1.0 - target.par_D(acoor));
                        }
                    }
                    acoor[1] =0;
                    for(acoor[0]=1 ;acoor[0]< target.par_RM.sizes[1]; acoor[0]++){
                        if (RM_deriv[acoor[0]][0] > RM_deriv[acoor[1]][0]) acoor[1] = acoor[0];
                    }
                    candidates[0].toZero();
                    candidates[0][acoor[1]] = 1.0;
                    if (candidates[0].deref_key(0) != curcol.deref_key(0)){
                        nbmax_diff[threadID][candidates[0].deref_key(0)]++;
                        nbmax_diff[threadID][curcol.deref_key(0)]--;
                        for(i = 0; i < data.data[coor[1]].getSize();i++){
                            acoor[0] = target.rowhid.data[ data.data[coor[1]].deref_key(i)].deref_key(0);
                            acoor[1] = curcol.deref_key(0);
                            nonzero_diff[threadID](acoor)--;
                            acoor[1] = candidates[0].deref_key(0);
                            nonzero_diff[threadID](acoor)++;
                        }
                    }
                    newhidden.memmoveColumn(coor[1], candidates[0]);*/
                }else{
                    checksafe = false;
                    if (curcol.getSize() == 0) {printf("Queried column %i is empty... impossible!\n",coor[1]); exit(1);}

                    if (coor[1] >= target.colhid.data.getSize()) exit(1);
                    // Step 1 cache line specific data
                    curproj[0] = target.par_RM * curcol;
                    curproj[3] = oldRM * curcol;
                    // Step 2 find direction, find 2 vector for R and M
                    RMS_deriv.toZero();




                    lls[2] = target.col_states_prior(curcol) * ((use_heavy_prior) ? data.data[coor[1]].getSize() : 1);
                    lls[3] = lls[2];
                    future_statistics[2].toZero(); future_statistic_trix[2].toZero();
                    d_scale[0].toZero();
                    d_scale[1].toZero();
                    dd_scale.toZero();
                    for(i=0;i<data.data[coor[1]].getSize();i++){
                        coor[0]= data.data[coor[1]].deref_key(i);
                        if (!ExOp::isValid(target.par_rowscale[coor[0]])) continue;
                        if (target.rowhid.getColumn(coor[0]).getSize() == 0) printf("dasize is empty!\n");
                        coeff = curproj[3].mkInnerProd(target.rowhid.getColumn(coor[0]));
                        coeff += target.par_rowscale[coor[0]];
                        coeff += target.par_colscale[coor[1]];
                        deflatedNBterm_wr_f_d2(ld.evilterms, data.data[coor[1]].deref(i)  , coeff[0], coeff[1]);
                        lls[3] += ld.evilterms[0];

                        coeff = curproj[0].mkInnerProd(target.rowhid.getColumn(coor[0]));
                        coeff += target.par_rowscale[coor[0]];
                        coeff += target.par_colscale[coor[1]];
                        deflatedNBterm_wr_f_d2(ld.evilterms, data.data[coor[1]].deref(i)  , coeff[0], coeff[1]);



                        dd_scale.data[0] += ld.evilterms[3]; dd_scale.data[1] += ld.evilterms[5]; dd_scale.data[2] += ld.evilterms[4];
                        d_scale[1][0] += ld.evilterms[1]; d_scale[1][1] += ld.evilterms[2];

                        lls[2] += ld.evilterms[0];

                        rever = target.par_RM.mkBackMult(target.rowhid.getColumn(coor[0]));

                        if (auto ite = rever.getIterator()) do{
                            RMS_deriv[ite()] += (*ite)[0] * ld.evilterms[1] + (*ite)[1] * ld.evilterms[2];
                        }while(ite++);

                        if (auto ite = target.rowhid.getColumn(coor[0]).getIterator()) do{

                            future_statistics[2][ite()][0] += ld.evilterm_alt.d[0] * (*ite);
                            future_statistics[2][ite()][1] += ld.evilterm_alt.d[1] * (*ite);
                            trix_input[0] = ld.evilterms[3] * (*ite);
                            trix_input[1] = ld.evilterms[4] * (*ite);
                            trix_input[2] = ld.evilterms[5] * (*ite) * 0.5;

                            if (auto ite2 = target.rowhid.getColumn(coor[0]).getIterator()) do{
                                if (ite2() > ite()) continue;
                                future_statistic_trix[2].data[ite2() + ((ite() * (ite() + 1)) >> 1)] += trix_input * (*ite2);
                            }while(ite2++);
                        }while(ite++);

                        rever_buf[i].toMemmove(rever);
                        dd_buf[i][0] = ld.evilterms[0];
                        dd_buf[i][1] = ld.evilterms[1];
                        dd_buf[i][2] = ld.evilterms[2];
                        dd_buf[i][3] = ld.evilterms[3];
                        dd_buf[i][4] = ld.evilterms[4];
                        dd_buf[i][5] = ld.evilterms[5];

                    } // dG/dz = c1f1 VM  dH/dz = c2 f2VR
                    //printf("cur coor %i/%i\n", coor[1], target.colhid.getNBcols());fflush(stdout);


                    prediction = 0.0f;
                    // Step 3 (line search) define set of candidates


                    sum[0] = 0.0;
                    prior_llterm[2] = 0.0;
                    candidates[0].toMemfree();

                    mergesuggestion = target.col_states_prior.findCandidateStep(candidates[0],curcol, wtorem, RMS_deriv.getSize(),RMS_deriv(),(use_heavy_prior) ? data.data[coor[1]].getSize() : 1);

                    // w has to be in 0-1 range
                    tosolve_num2[0] = 0; //(target.colhid_prior.mkInnerMult(candidates[0],curcol) + target.colhid_prior_constant) * ( 2.0 * ((tmptmp) ? data.data[coor[1]].getSize() : 1)/ prior_llterm[0]);
                    tosolve_num2[1] = d_scale[1][0];
                    tosolve_num2[2] = d_scale[1][1];
                    tosolve_den2.data[1] = 0.0; tosolve_den2.data[3] = 0.0;
                    tosolve_den2.data[2] = dd_scale.data[0];
                    tosolve_den2.data[5] = dd_scale.data[2];
                    tosolve_den2.data[4] = dd_scale.data[1];
                    tosolve_den2.data[0] = 0.0; //(target.colhid_prior.mkInnerMult(candidates[0])+ target.colhid_prior_constant) * ( 2.0 * ((tmptmp) ? data.data[coor[1]].getSize() : 1) / prior_llterm[0]);
                    for(i=0;i<data.data[coor[1]].getSize();i++){
                        coor[0]= data.data[coor[1]].deref_key(i);
                        if (!ExOp::isValid(target.par_rowscale[coor[0]])) continue;
                        coeff = rever_buf[i].mkInnerProd(candidates[0]);
                        tosolve_num2[0] += coeff[0] * dd_buf[i][1]+ coeff[1] * dd_buf[i][2];
                        tosolve_den2.data[0] += dd_buf[i][3] * coeff[0] * coeff[0] + dd_buf[i][4] * coeff[1] * coeff[1] + dd_buf[i][5] * coeff[0] * coeff[1];
                        tosolve_den2.data[1] += dd_buf[i][3] * coeff[0] + dd_buf[i][5] * coeff[1] * 0.5;
                        tosolve_den2.data[3] += dd_buf[i][4] * coeff[1] + dd_buf[i][5] * coeff[0] * 0.5;
                    }

                    tosolve_den2_tmp = tosolve_den2;
                    tosolve_den2.toEigenTransform(funabsoft);
                    solution = tosolve_den2 * tosolve_num2;

                    expexp = solution[0] * tosolve_num2[0] + solution[1] * tosolve_num2[1] + solution[2] * tosolve_num2[2];
                    expexp2 = tosolve_den2_tmp.data[0] * solution[0] * solution[0] * 0.5;
                    expexp2 += tosolve_den2_tmp.data[2] * solution[1] * solution[1] * 0.5;
                    expexp2 += tosolve_den2_tmp.data[5] * solution[2] * solution[2] * 0.5;
                    expexp2 += tosolve_den2_tmp.data[1] * solution[0] * solution[1];
                    expexp2 += tosolve_den2_tmp.data[3] * solution[0] * solution[2];
                    expexp2 += tosolve_den2_tmp.data[4] * solution[1] * solution[2];


                    dd_scale.toEigenTransform(funabsoft);
                    if ((solution[0] <= 0.0)||(solution[0] >= 1.0)) { // solution is outside allowed range... try dropping dimension or just update scale only
                        candidates[1] = curcol;
                        if (mergesuggestion[0] != 0xFFFFFFFF) {
                            tosolve_den2.data[1] = tosolve_den2.data[3] = 0.0f;
                            for(i=0;i<data.data[coor[1]].getSize();i++){
                                coor[0] = data.data[coor[1]].deref_key(i);
                                if (!ExOp::isValid(target.par_rowscale[coor[0]])) continue;
                                coeff = rever_buf[i][mergesuggestion[1]] - rever_buf[i][mergesuggestion[0]];
                                tosolve_den2.data[1] += dd_buf[i][3] * coeff[0] + dd_buf[i][5] * coeff[1] * 0.5;
                                tosolve_den2.data[3] += dd_buf[i][4] * coeff[1] + dd_buf[i][5] * coeff[0] * 0.5;
                            }
                            d_scale[1][0] = tosolve_num2[1] + tosolve_den2.data[1] * curcol[mergesuggestion[0]];
                            d_scale[1][1] = tosolve_num2[2] + tosolve_den2.data[3] * curcol[mergesuggestion[0]];
                            candidates[1][mergesuggestion[1]] += curcol[mergesuggestion[0]];
                            candidates[1].erase(mergesuggestion[0]);
                        }
                        d_scale[1] = dd_scale * d_scale[1];
                        candidates[0] += curcol;
                        candidates[0].erase(wtorem);
                        d_scale[0][0] = tosolve_num2[1] + tosolve_den2.data[1]; // deritive evaluated at pt
                        d_scale[0][1] = tosolve_num2[2] + tosolve_den2.data[3]; // deritive evaluated at pt
                    }else{ // all set!
                        d_scale[1][0] = solution[1];
                        d_scale[1][1] = solution[2];
                        candidates[1] = curcol + (candidates[0] * solution[0]);
                        candidates[1].keysort();
                        sum[0] = 0.0f;
                        if (auto ite = candidates[1].getIterator()) do{
                            sum[0] += *ite;
                        }while(ite++);
                        if (fabs(sum[0] - 1.0) > 0.1) {printf("hidden state norm is hh %e! fact %e \n", sum[0], solution[0]);
                            candidates[0].show();
                            candidates[1].show();
                            curcol.show();
                            exit(1);
                        }
                        if ((mergesuggestion[0] == 0xFFFFFFFF)||((rand() & 128) != 0)) {
                            candidates[0] += curcol;
                            candidates[0].erase(wtorem);
                            d_scale[0][0] = tosolve_num2[1] + tosolve_den2.data[1]; // deritive evaluated at pt
                            d_scale[0][1] = tosolve_num2[2] + tosolve_den2.data[3]; // deritive evaluated at pt
                        }else{
                            tosolve_den2.data[1] = tosolve_den2.data[3] = 0.0f;
                            for(i=0;i<data.data[coor[1]].getSize();i++){
                                coor[0] = data.data[coor[1]].deref_key(i);
                                if (!ExOp::isValid(target.par_rowscale[coor[0]])) continue;
                                coeff = rever_buf[i][mergesuggestion[1]] - rever_buf[i][mergesuggestion[0]];
                                tosolve_den2.data[1] += dd_buf[i][3] * coeff[0] + dd_buf[i][5] * coeff[1] * 0.5;
                                tosolve_den2.data[3] += dd_buf[i][4] * coeff[1] + dd_buf[i][5] * coeff[0] * 0.5;
                            }
                            d_scale[0][0] = tosolve_num2[1] + tosolve_den2.data[1] * curcol[mergesuggestion[0]];
                            d_scale[0][1] = tosolve_num2[2] + tosolve_den2.data[3] * curcol[mergesuggestion[0]];
                            candidates[0] = curcol;
                            candidates[0][mergesuggestion[1]] += curcol[mergesuggestion[0]];
                            candidates[0].erase(mergesuggestion[0]);
                        }
                    }
                    d_scale[0] = dd_scale * d_scale[0];
                    candidates[0].keysort();
                    sum[0] = 0.0f;
                    if (auto ite = candidates[0].getIterator()) {
                        do{
                            sum[0] += *ite;
                        }while(ite++);
                        if (fabs(sum[0] - 1.0) > 0.1) {printf("hidden state norm is hh %e!\n", sum[0]);
                            candidates[0].show();
                            curcol.show();
                        exit(1);}
                    }else {printf("got empty candidate!\n"); exit(1);}


                    // find intersection on line
                    // new scale!
                    //printf("cur coor %i/%i\n", coor[1], target.colhid.getNBcols());fflush(stdout);


                    // Step 4 Evaluate Quantidates

                    lld[2] = 0.0f;
                    tosolve_num2[0] = 1.0 / curcol.getSize();
                    if (auto ite = curcol.getIterator()) do{
                        acoor[1] = ite();
                        for(acoor[0] =0;acoor[0]< target.par_RM.sizes[0]; acoor[0]++){
                            coor[0] = acoor[0];
                            if (freq_rowh[coor[0]] <= col_freq_hrow(coor)- 0.001) thrbase.msgs.insert(to_string(coor[0]) + string("got negative!") + to_string((double)col_freq_hrow(coor)));
                            lld[2] += tosolve_num2[0] * (col_freq_hrow(coor) * log(target.par_D(acoor)) + (freq_rowh[coor[0]] - col_freq_hrow(coor)) * log(1.0 - target.par_D(acoor)));
                        }
                    }while(ite++);
                    lld[2] *= dropout_likelihood_factor;
                    if (!ExOp::isValid(lld[2])){
                        thrbase.msgs.insert(string("Default Dropout likelihood is NA"));
                        thrbase.msgs.insert(string("ite current size is") + to_string(curcol.getSize()));
                    }

                    insist =0;
                    maindoutcmp = false;
                    while(true){

                        for(j=0;j<2;j++){
                            future_statistics[j].toZero(); future_statistic_trix[j].toZero();
                            curproj[j] = target.par_RM * candidates[j];
                            lls[j] = target.col_states_prior(candidates[j]) * ((use_heavy_prior) ? data.data[coor[1]].getSize() : 1);
                            if (candidates[j].hasSameSparsity(curcol)) lld[j] = lld[2];
                            else{
                                lld[j] = 0.0f;
                                tosolve_num2[0] = 1.0 / candidates[j].getSize();
                                if (auto ite = candidates[j].getIterator()) do{
                                    acoor[1] = ite();
                                    for(acoor[0] =0;acoor[0]< target.par_RM.sizes[0]; acoor[0]++){
                                        coor[0] = acoor[0];
                                        lld[j] += tosolve_num2[0] * (col_freq_hrow(coor) * log(target.par_D(acoor)) + (freq_rowh[coor[0]] - col_freq_hrow(coor)) * log(1.0 - target.par_D(acoor)));
                                    }
                                }while(ite++);
                                if (!ExOp::isValid(lld[j])) thrbase.msgs.insert(string("Canditate dropoutLL is invalid! ")+ to_string(j));
                                lld[j] *= dropout_likelihood_factor;
                            }
                            if ((candidates[1].getSize() == 0)||(insist > 0)) break;
                        }
                        for(i=0;i<data.data[coor[1]].getSize();i++){
                            coor[0]= data.data[coor[1]].deref_key(i);
                            if (!ExOp::isValid(target.par_rowscale[coor[0]])) continue;
                            for(j=0;j<2;j++){
                                // first check if scale update is dangerous
                                coeff = target.par_colscale[coor[1]];
                                coeff += d_scale[j];
                                if ((tmp = (coeff[0]*coeff[0] + coeff[1]*coeff[1])) > 100.0){ // too far, send to boundary
                                    tmp = pow(tmp, -0.5);
                                    d_scale[j][0] = coeff[0] * tmp;
                                    d_scale[j][1] = coeff[1] * tmp;
                                    d_scale[j] -= target.par_colscale[coor[1]];
                                }

                                coeff = curproj[j].mkInnerProd(target.rowhid.getColumn(coor[0]));
                                coeff += target.par_rowscale[coor[0]];
                                coeff += target.par_colscale[coor[1]];
                                coeff += d_scale[j];


                                deflatedNBterm_wr_f_d2(ld.evilterms, data.data[coor[1]].deref(i), coeff[0], coeff[1]);
                                if (auto ite = target.rowhid.getColumn(coor[0]).getIterator()) do{
                                    future_statistics[j][ite()][0] += ld.evilterm_alt.d[0] * (*ite);
                                    future_statistics[j][ite()][1] += ld.evilterm_alt.d[1] * (*ite);
                                    trix_input[0] = ld.evilterms[3] * (*ite);
                                    trix_input[1] = ld.evilterms[4] * (*ite);
                                    trix_input[2] = ld.evilterms[5] * (*ite) * 0.5;
                                    if (auto ite2 = target.rowhid.getColumn(coor[0]).getIterator()) do{
                                        if (ite2() > ite()) continue;
                                        future_statistic_trix[j].data[ite2() + ((ite() * (ite() + 1)) >> 1)] += trix_input * (*ite2);
                                    }while(ite2++);
                                }while(ite++);
                                lls[j] += ld.evilterms[0];
                            }


                        } // H(z) = fcfc vMMv





                        /*if (checksafe){
                            double diffdiff = (log(target.colhid_prior.mkInnerMult(candidates[0])+ target.colhid_prior_constant) - log(target.colhid_prior.mkInnerMult(curcol)+ target.colhid_prior_constant)) * ((tmptmp) ? data.data[coor[1]].getSize() : 1);
                            printf("%e vs %e predicted\n", lls[0] - lls[2], expexp + expexp2);
                            printf("quad term: %e vs %e\n", lls[0] - lls[2] - expexp, expexp2);
                        }*/

                        // Step 4.1 Append LL from dropouts  TODO, ignore drop out for now
                        /*if (haschange){
                            for(i=0;i< nb_maxrowh.getSize();i++{


                            }
                        }*/

                        // Step 5 select best or do not update
                        //tmp = d_scale_tmp[0] * d_scale[0] + d_scale_tmp[1] * d_scale[1];
                        //tmp += (dd_scale_tmp.data[0] * d_scale[0] * d_scale[0] + dd_scale_tmp.data[2] * d_scale[1] * d_scale[1])* 0.5 + dd_scale_tmp.data[1] * d_scale[0] * d_scale[1];
                        //printf("%e incr (vs%e)\n", lls[0] - lls[3], tmp);

                    //    for(coor[0]=0;coor[0] < new_nbmax_h.getSize();coor[0]++) new_nbmax_h(coor) = row_nbmax_hcol(coor);


                        j = ((candidates[1].getSize() == 0)||(lls[0] + lld[0] > lls[1] + lld[1])) ? 0 : 1;
                        j = (lls[j] + lld[j] > lls[2] + lld[2]) ? j : 2;
                        if (j != 2){
                            if (candidates[j].hasSameSparsity(curcol)) update_type[coor[1]] = 'S';
                            else {
                                update_type[coor[1]] = 'D';
                                // registering change...

                                tosolve_num2[0] = 1.0 / curcol.getSize();
                                tosolve_num2[1] = 1.0 / candidates[j].getSize();
                                if (auto ite = curcol.getIterator()) do{freq_diff[threadID][ite()] -= tosolve_num2[0];} while(ite++);
                                if (auto ite = candidates[j].getIterator()) do{freq_diff[threadID][ite()] += tosolve_num2[1];} while(ite++);

                                for(i = 0; i < data.data[coor[1]].getSize();i++){
                                    tosolve_num2[0] = 1.0 / (curcol.getSize() *  target.rowhid.data[ data.data[coor[1]].deref_key(i)].getSize());
                                    tosolve_num2[1] = 1.0 / (candidates[j].getSize() *  target.rowhid.data[ data.data[coor[1]].deref_key(i)].getSize());
                                    if (auto ite2 = target.rowhid.data[ data.data[coor[1]].deref_key(i)].getIterator()) do{
                                        acoor[0] = ite2();
                                        if (auto ite = curcol.getIterator()) do{acoor[1] = ite(); nonzero_freq_diff[threadID](acoor) -= tosolve_num2[0];} while(ite++);
                                        if (auto ite = candidates[j].getIterator()) do{acoor[1] = ite(); nonzero_freq_diff[threadID](acoor) += tosolve_num2[1];} while(ite++);
                                    }while(ite2++);
                                }
                                nb_flags[threadID][2]++;
                                /*
                                nbmax_diff[threadID][candidates[j].deref_key(0)]++;
                                nbmax_diff[threadID][curcol.deref_key(0)]--;
                                for(i = 0; i < data.data[coor[1]].getSize();i++){
                                    acoor[0] = target.rowhid.data[ data.data[coor[1]].deref_key(i)].deref_key(0);
                                    acoor[1] = curcol.deref_key(0);
                                    nonzero_diff[threadID](acoor)--;
                                    acoor[1] = candidates[j].deref_key(0);
                                    nonzero_diff[threadID](acoor)++;
                                }*/
                            }
                        }

                        if (j != 2) break;
                        if (insist == 0) {
                            if (candidates[1] == curcol) {update_type[coor[1]] = 'F'; thrbase.msgs.insert(string("Failed to update scale alone!\n")); break;}
                            candidates[0].toMemmove(candidates[1]);
                        }
                        if (insist++ > 10) {update_type[coor[1]] = 'F'; break;}

                        if (auto ite = curcol.getIterator()) do{
                            candidates[0][ite()] += *ite;
                        }while(ite++);

                        if (auto ite = candidates[0].getIterator()) do{
                            *ite *= 0.5;
                        }while(ite++);
                        d_scale[0] *= 0.5;
                    } // insist loop end
                    nb_flags[threadID][1] += insist;

                    // Step 6 save scope for parameter update
                    /*ll[threadID][0] -= target.col_states_prior(curcol);
                    ll[threadID][1] -= target.col_states_prior(candidates[j]);
                    ll[threadID][4] += target.col_states_prior(curcol);
                    ll[threadID][5] += target.col_states_prior(candidates[j]);*/
                    newscale[coor[1]] = target.par_colscale[coor[1]];
                    if (j == 2) { // failed to improve...
                        //printf("col %i failed  %e %e\n", coor[1], target.par_colscale[coor[1]][0], target.par_colscale[coor[1]][1] );
                        candidates[0] = curcol;
                        newhidden.memmoveColumn(coor[1], candidates[0]);
                        nb_flags[threadID][0]++;
                    }else{
                        newhidden.memmoveColumn(coor[1], candidates[j]);
                        newscale[coor[1]] += d_scale[j];
                    }

                    ll[threadID][0] += lls[2];
                    ll[threadID][1] += lls[j];
                    ll[threadID][2] += lld[2];
                    ll[threadID][3] += lld[j];
                    ll[threadID][4] += lls[3];

                    /*if ((coor[1] % 100) == 57){
                        thrbase.msgs.insert( to_string((int) coor[1]) + string("LLs[j]") + to_string((double)lls[j]));
                        thrbase.msgs.insert( to_string((int) coor[1]) + string("LLs[3]") + to_string((double)lls[3]));
                    }*/
                    candidates[1].toMemfree();

                    if (auto ite =newhidden.data[coor[1]].getIterator()) do {if (ite() >= target.par_RM.sizes[1]) {printf("IMPOfsdSSdddIBLBLBLB!!!\n"); curcol.show(); exit(1);} }while(ite++);



                    future_factor.toMemfree();
                    if (auto ite = newhidden.data[coor[1]].getIterator()) do{
                        future_factor[ite()][0] = *ite; future_factor[ite()][1] = *ite;
                    }while(ite++);

                    dll_dMR2[threadID].toAddOuterProd(future_statistics[j], future_factor);

                    uint32_t dashift = target.par_RM.sizes[0] * target.par_RM.sizes[1];
                    if (auto ite = newhidden.data[coor[1]].getIterator()) do{
                        if (auto ite2 = newhidden.data[coor[1]].getIterator()) do{
                            if (ite2() > ite()) continue;
                            tmp = (*ite2) * (*ite);
                            uint32_t k;
                            for(i = 0,k=0; i< target.par_RM.sizes[0];i++){
                                uint32_t l = i + ite() * target.par_RM.sizes[0];

                                double* currow_A = dd_massive[threadID].data + (((l * (l +1)) >> 1) + ite2() * target.par_RM.sizes[0]);
                                double* currow_B = dd_massive[threadID].data + ((((l + dashift) * (l + dashift +1)) >> 1) + ite2() * target.par_RM.sizes[0]);
                                for(l=0;l<=i;l++,k++){
                                    currow_A[l] += future_statistic_trix[j].data[k][0] * tmp;
                                    currow_B[l] += future_statistic_trix[j].data[k][2] * tmp;
                                    currow_B[l + dashift] += future_statistic_trix[j].data[k][1] * tmp;
                                }
                            }
                            //printf("%i is k  %i\n", k,  future_statistic_trix[j].getSize());
                        }while(ite2++);
                    }while(ite++);
                }
            }
            lls[0] =0;
            for(coor[1]=rowrange[threadID];coor[1] < rowrange[threadID+1];coor[1]++){
                if (trdata.getColSize(coor[1]) != 0) lls[0] += target.row_states_prior(target.rowhid.getColumn(coor[1])) * ((use_heavy_prior) ? trdata.data[coor[1]].getSize() : 1);
            }
            if (!ExOp::isValid(lls[0])) thrbase.msgs.insert( string("thread")+  to_string(threadID) +  string("got invalid likelihood"));
            ll[threadID][0] += lls[0];
            ll[threadID][1] += lls[0];
            ll[threadID][4] += lls[0];
        break;
        case 1: // time for rows
            RMS_deriv.setSize(target.par_RM.sizes[0]);
            rever_buf.setSize(target.colhid.getNBcols());
            dd_buf.setSize(target.colhid.getNBcols());

            for(i=0;i<3;i++) {
                future_statistics[i].setSize(target.par_RM.sizes[1]);
                future_statistic_trix[i].setSize(target.par_RM.sizes[1]);
            }
            ExOp::toZero(sum);

            for(coor[1]=rowrange[threadID];coor[1] < rowrange[threadID+1];coor[1]++){ thrbase.updateProgress(threadID);
                //printf("d%ie%if\n", threadID, coor[1]);
                if (trdata.data[coor[1]].getSize() > target.colhid.getNBcols()) {printf("Holly Molly! %i and %i",trdata.data[coor[1]].getSize(), target.colhid.getNBcols());exit(1);}
                //printf("cur coor %i/%i\n", coor[1], target.colhid.getNBcols());fflush(stdout);
                const SparseTuple<double>& curcol = target.rowhid.getColumn(coor[1]);
                if (!ExOp::isValid(target.par_rowscale[coor[1]][0])){
                    newscale[coor[1]] = target.par_rowscale[coor[1]];

                    if (trdata.getColSize(coor[1]) == 0) {
                            update_type[coor[1]] = 'S'; candidates[0] = curcol;newhidden.memmoveColumn(coor[1], candidates[0]); continue;
                    } // nothing to see...
                    // column only has 1s, cannot be used to update R or M, drop out only, pure state guarrantied
                    if (curcol.getSize() != 1) {
                        printf("half-unsuable row%i (with %i non-zero counts)should have 1 dim only... got %i\n",coor[1], trdata.getColSize(coor[1]), curcol.getSize());
                        for(i=0;i<100;i++) printf("%i\t%i\n", trdata.getColSize(coor[1]), target.rowhid.getColumn(coor[1]).getSize());
                    exit(1);}

                    // try 2 random alternatives
                     candidates[0].toMemfree();
                     candidates[1].toMemfree();
                     i = rand() % (target.par_RM.sizes[0] - 1);
                     if (i >= curcol.deref_key(0)) i++;
                     candidates[0][i] = 1.0;
                     lls[0] = target.row_states_prior(candidates[0]) * ((use_heavy_prior) ? trdata.data[coor[1]].getSize() : 1);
                     lls[2] = target.row_states_prior(curcol) * ((use_heavy_prior) ? trdata.data[coor[1]].getSize() : 1);
                     if (target.par_RM.sizes[0] > 2){
                         j = rand() % (target.par_RM.sizes[0] - 2);
                         if (j >= curcol.deref_key(0)) j++;
                         if (j >= i) j++;
                        candidates[1][j] = 1.0;
                        lls[1] = target.row_states_prior(candidates[1]) * ((use_heavy_prior) ? trdata.data[coor[1]].getSize() : 1);
                     }else j=i;
                     ExOp::toZero(lld);

                     acoor[0] = curcol.deref_key(0);
                     for(acoor[1]=0;acoor[1]< target.par_RM.sizes[1]; acoor[1]++){
                         coor[0] = acoor[1];
                         if (freq_colh[coor[0]] <= row_freq_hcol(coor)- 0.001) thrbase.msgs.insert(string("got negative!") + to_string(freq_colh[coor[0]] - row_freq_hcol(coor)));
                         lld[2] += row_freq_hcol(coor) * log(target.par_D(acoor)) + (freq_colh[coor[0]] - row_freq_hcol(coor)) * log(1.0 - target.par_D(acoor));
                     }
                     acoor[0] = candidates[0].deref_key(0);
                     for(acoor[1]=0;acoor[1]< target.par_RM.sizes[1]; acoor[1]++){
                         coor[0] = acoor[1];
                         if (freq_colh[coor[0]] <= row_freq_hcol(coor)- 0.001) thrbase.msgs.insert(string("got negative!") + to_string(freq_colh[coor[0]] - row_freq_hcol(coor)));
                         lld[0] += row_freq_hcol(coor) * log(target.par_D(acoor)) + (freq_colh[coor[0]] - row_freq_hcol(coor)) * log(1.0 - target.par_D(acoor));
                     }
                     acoor[0] = j;
                     for(acoor[1]=0;acoor[1]< target.par_RM.sizes[1]; acoor[1]++){
                         coor[0] = acoor[1];
                         if (freq_colh[coor[0]] <= row_freq_hcol(coor)- 0.001) thrbase.msgs.insert(string("got negative!") + to_string(freq_colh[coor[0]] - row_freq_hcol(coor)));
                         lld[1] += row_freq_hcol(coor) * log(target.par_D(acoor)) + (freq_colh[coor[0]] - row_freq_hcol(coor)) * log(1.0 - target.par_D(acoor));
                     }
                     lld[0] *= dropout_likelihood_factor;
                     lld[1] *= dropout_likelihood_factor;
                     lld[2] *= dropout_likelihood_factor;
                     j = ((candidates[1].getSize() == 0)||(lls[0] + lld[0] > lls[1] + lld[1])) ? 0 : 1;
                     j = (lls[j] + lld[j] > lls[2] + lld[2]) ? j : 2;
                     if (j == 2){
                        candidates[0] = curcol; newhidden.memmoveColumn(coor[1], candidates[0]);
                        update_type[coor[1]] = 'S';
                     }else{
                        newhidden.memmoveColumn(coor[1], candidates[j]);
                        update_type[coor[1]] = 'D';

                        freq_diff[threadID][curcol.deref_key(0)] -= 1.0;
                        freq_diff[threadID][candidates[j].deref_key(0)] += 1.0;
                        for(i = 0; i < trdata.data[coor[1]].getSize();i++){
                            tosolve_num2[0] = 1.0 / (target.colhid.data[trdata.data[coor[1]].deref_key(i)].getSize());
                            if (auto ite2 = target.colhid.data[trdata.data[coor[1]].deref_key(i)].getIterator()) do{
                                acoor[1] = ite2();
                                acoor[0] = curcol.deref_key(0);        nonzero_freq_diff[threadID](acoor) -= tosolve_num2[0];
                                acoor[0] = candidates[j].deref_key(0); nonzero_freq_diff[threadID](acoor) += tosolve_num2[0];
                            }while(ite2++);
                        }
                        nb_flags[threadID][2]++;
                     }
                    ll[threadID][0] += lls[2];
                    ll[threadID][1] += lls[j];
                    ll[threadID][2] += lld[2];
                    ll[threadID][3] += lld[j];
                    ll[threadID][4] += lls[2];

                    candidates[1].toMemfree();
                    candidates[0].toMemfree();
                    continue;
                    // put in stone... for now
                    /*
                    RM_deriv.toZero();

                    for(coor[0]=0;coor[0]<row_nbmax_hcol.dims[0];coor[0]++){acoor[1] = coor[0];
                        for(acoor[0] =0;acoor[0]< target.par_RM.sizes[0]; acoor[0]++){
                            RM_deriv[acoor[0]][0] += row_nbmax_hcol(coor) * log(target.par_D(acoor));
                            RM_deriv[acoor[0]][0] += (nb_maxcolh[coor[0]] - row_nbmax_hcol(coor)) * log(1.0 - target.par_D(acoor));
                        }
                    }
                    acoor[1] =0;
                    for(acoor[0]=0 ;acoor[0]< target.par_RM.sizes[0]; acoor[0]++){
                        if (RM_deriv[acoor[0]][0] > RM_deriv[acoor[1]][0]) acoor[1] = acoor[0];
                    }
                    candidates[0].toZero();
                    candidates[0][acoor[1]] = 1.0;

                    if (candidates[0].deref_key(0) != curcol.deref_key(0)){
                        nbmax_diff[threadID][candidates[0].deref_key(0)]++;
                        nbmax_diff[threadID][curcol.deref_key(0)]--;
                        for(i = 0; i < trdata.data[coor[1]].getSize();i++){
                            acoor[1] = target.colhid.data[trdata.data[coor[1]].deref_key(i)].deref_key(0);
                            acoor[0] = curcol.deref_key(0);
                            nonzero_diff[threadID](acoor)--;
                            acoor[0] = candidates[0].deref_key(0);
                            nonzero_diff[threadID](acoor)++;
                        }
                    }
                    newhidden.memmoveColumn(coor[1], candidates[0]);*/
                }else{
                    if (curcol.getSize() == 0) {printf("Queried row %i is empty... impossible!\n",coor[1]); exit(1);}

                    // Step 1 cache line specific data
                    curproj[0] = target.par_RM.mkBackMult(curcol);
                    curproj[3] = oldRM.mkBackMult(curcol);
                    //dll_dH.toZero();
                    //dll_dC.toZero();

                    // Step 2 find direction, find 2 vector for R and M
                    RMS_deriv.toZero();
                    lls[2] = target.row_states_prior(curcol) * ((use_heavy_prior) ? trdata.data[coor[1]].getSize() : 1);
                    lls[3] = lls[2];
                    future_statistics[2].toZero();future_statistic_trix[2].toZero();
                    d_scale[0].toZero();
                    d_scale[1].toZero();

                    dd_scale.toZero();
                   // printf("loopya\n"); fflush(stdout);
                    for(i=0;i<trdata.data[coor[1]].getSize();i++){
                        coor[0]= trdata.data[coor[1]].deref_key(i);
                        if (!ExOp::isValid(target.par_colscale[coor[0]][0])) continue; // can do better...
                        coeff = curproj[3].mkInnerProd(target.colhid.getColumn(coor[0]) );
                        coeff += target.par_colscale[coor[0]];
                        coeff += target.par_rowscale[coor[1]];
                        deflatedNBterm_wr_f_d2(ld.evilterms, trdata.data[coor[1]].deref(i), coeff[0], coeff[1]);
                        lls[3] += ld.evilterms[0];

                        coeff = curproj[0].mkInnerProd(target.colhid.getColumn(coor[0]) );
                        coeff += target.par_colscale[coor[0]];
                        coeff += target.par_rowscale[coor[1]];
                        deflatedNBterm_wr_f_d2(ld.evilterms, trdata.data[coor[1]].deref(i), coeff[0], coeff[1]);
                        dd_scale.data[0] += ld.evilterms[3]; dd_scale.data[1] += ld.evilterms[5]; dd_scale.data[2] += ld.evilterms[4];
                        d_scale[1][0] += ld.evilterms[1]; d_scale[1][1] += ld.evilterms[2];
                        lls[2] += ld.evilterms[0];

                        rever = target.par_RM * target.colhid.getColumn(coor[0]);
                        if (auto ite = rever.getIterator()) do{
                            RMS_deriv[ite()] += (*ite)[0] * ld.evilterms[1] + (*ite)[1] * ld.evilterms[2];
                        }while(ite++);


                        if (auto ite = target.colhid.getColumn(coor[0]).getIterator()) do{
                            future_statistics[2][ite()][0] += ld.evilterm_alt.d[0] * (*ite);
                            future_statistics[2][ite()][1] += ld.evilterm_alt.d[1] * (*ite);
                            trix_input[0] = ld.evilterms[3] * (*ite);
                            trix_input[1] = ld.evilterms[4] * (*ite);
                            trix_input[2] = ld.evilterms[5] * (*ite) * 0.5;
                            if (auto ite2 = target.colhid.getColumn(coor[0]).getIterator()) do{
                                if (ite2() > ite()) continue;
                                future_statistic_trix[2].data[ite2() + ((ite() * (ite() + 1)) >> 1)] += trix_input * (*ite2);
                            }while(ite2++);

                        }while(ite++);

                        rever_buf[i].toMemmove(rever);
                        dd_buf[i][0] = ld.evilterms[0];
                        dd_buf[i][1] = ld.evilterms[1];
                        dd_buf[i][2] = ld.evilterms[2];
                        dd_buf[i][3] = ld.evilterms[3];
                        dd_buf[i][4] = ld.evilterms[4];
                        dd_buf[i][5] = ld.evilterms[5];
                    }
                    // Step 3 (line search) define set of candidates
                    sum[0] = 0.0; prior_llterm[2] = 0.0;
                    candidates[0].toMemfree();
                    mergesuggestion = target.row_states_prior.findCandidateStep(candidates[0],curcol, wtorem, RMS_deriv.getSize(),RMS_deriv(),(use_heavy_prior) ? data.data[coor[1]].getSize() : 1);
                    coor[0] =0;
                    tosolve_num2[0] = 0.0;//(target.rowhid_prior.mkInnerMult(candidates[0],curcol)+ target.rowhid_prior_constant) * ( 2.0 * ((tmptmp) ? trdata.data[coor[1]].getSize() : 1)  / prior_llterm[0]);
                    tosolve_num2[1] = d_scale[1][0];
                    tosolve_num2[2] = d_scale[1][1];
                    tosolve_den2.data[1] = 0.0; tosolve_den2.data[3] = 0.0;
                    tosolve_den2.data[2] = dd_scale.data[0];
                    tosolve_den2.data[5] = dd_scale.data[2];
                    tosolve_den2.data[4] = dd_scale.data[1];
                    tosolve_den2.data[0] = 0.0; // (target.rowhid_prior.mkInnerMult(candidates[0])+ target.rowhid_prior_constant) * ( 2.0 * ((tmptmp) ? trdata.data[coor[1]].getSize() : 1) / prior_llterm[0]);

                    for(i=0;i<trdata.data[coor[1]].getSize();i++){
                        coor[0]= trdata.data[coor[1]].deref_key(i);
                        if (!ExOp::isValid(target.par_colscale[coor[0]])) continue;
                        coeff = rever_buf[i].mkInnerProd(candidates[0]);
                        tosolve_num2[0] += coeff[0] * dd_buf[i][1] + coeff[1] * dd_buf[i][2];
                        tosolve_den2.data[0] += dd_buf[i][3] * coeff[0] * coeff[0] + dd_buf[i][4] * coeff[1] * coeff[1] + dd_buf[i][5] * coeff[0] * coeff[1];
                        tosolve_den2.data[1] += dd_buf[i][3] * coeff[0] + dd_buf[i][5] * coeff[1] * 0.5;
                        tosolve_den2.data[3] += dd_buf[i][4] * coeff[1] + dd_buf[i][5] * coeff[0] * 0.5;
                    }

                    tosolve_den2_tmp = tosolve_den2;
                    tosolve_den2.toEigenTransform(funabsoft);
                    solution = tosolve_den2 * tosolve_num2;

                    expexp = solution[0] * tosolve_num2[0] + solution[1] * tosolve_num2[1] + solution[2] * tosolve_num2[2];
                    expexp2 = tosolve_den2_tmp.data[0] * solution[0] * solution[0] * 0.5;
                    expexp2 += tosolve_den2_tmp.data[2] * solution[1] * solution[1] * 0.5;
                    expexp2 += tosolve_den2_tmp.data[5] * solution[2] * solution[2] * 0.5;
                    expexp2 += tosolve_den2_tmp.data[1] * solution[0] * solution[1];
                    expexp2 += tosolve_den2_tmp.data[3] * solution[0] * solution[2];
                    expexp2 += tosolve_den2_tmp.data[4] * solution[1] * solution[2];


                    dd_scale.toEigenTransform(funabsoft);

                    if ((solution[0] <= 0.0)||(solution[0] >= 1.0)) { // solution is outside allowed range... try to update scale only
                        candidates[1] = curcol;
                        if (mergesuggestion[0] != 0xFFFFFFFF) {
                            tosolve_den2.data[1] = tosolve_den2.data[3] = 0.0f;
                            for(i=0;i<trdata.data[coor[1]].getSize();i++){
                                coor[0]= trdata.data[coor[1]].deref_key(i);
                                if (!ExOp::isValid(target.par_colscale[coor[0]])) continue;
                                coeff = rever_buf[i][mergesuggestion[1]] - rever_buf[i][mergesuggestion[0]];
                                tosolve_den2.data[1] += dd_buf[i][3] * coeff[0] + dd_buf[i][5] * coeff[1] * 0.5;
                                tosolve_den2.data[3] += dd_buf[i][4] * coeff[1] + dd_buf[i][5] * coeff[0] * 0.5;
                            }
                            d_scale[1][0] = tosolve_num2[1] + tosolve_den2.data[1] * curcol[mergesuggestion[0]];
                            d_scale[1][1] = tosolve_num2[2] + tosolve_den2.data[3] * curcol[mergesuggestion[0]];
                            candidates[1][mergesuggestion[1]] += curcol[mergesuggestion[0]];
                            candidates[1].erase(mergesuggestion[0]);
                        }
                        d_scale[1] = dd_scale * d_scale[1];
                        candidates[0] += curcol;
                        candidates[0].erase(wtorem);
                        d_scale[0][0] = tosolve_num2[1] + tosolve_den2.data[1];
                        d_scale[0][1] = tosolve_num2[2] + tosolve_den2.data[3];
                    }else{ // all set!
                        d_scale[1][0] = solution[1];
                        d_scale[1][1] = solution[2];
                        candidates[1] = curcol + (candidates[0] * solution[0]);
                        candidates[1].keysort();
                        sum[0] = 0.0f;
                        if (auto ite = candidates[1].getIterator()) do{
                            sum[0] += *ite;
                        }while(ite++);
                        if (fabs(sum[0] - 1.0) > 0.1) {printf("hidden state norm is hh %e! fact %e \n", sum[0], solution[0]);
                            candidates[0].show();
                            candidates[1].show();
                            curcol.show();
                            exit(1);
                        }
                        if ((mergesuggestion[0] == 0xFFFFFFFF)||((rand() & 128) != 0)) {
                            candidates[0] += curcol;
                            candidates[0].erase(wtorem);
                            d_scale[0][0] = tosolve_num2[1] + tosolve_den2.data[1];
                            d_scale[0][1] = tosolve_num2[2] + tosolve_den2.data[3];
                        }else{
                            tosolve_den2.data[1] = tosolve_den2.data[3] = 0.0f;
                            for(i=0;i<trdata.data[coor[1]].getSize();i++){
                                coor[0]= trdata.data[coor[1]].deref_key(i);
                                if (!ExOp::isValid(target.par_colscale[coor[0]])) continue;
                                coeff = rever_buf[i][mergesuggestion[1]] - rever_buf[i][mergesuggestion[0]];
                                tosolve_den2.data[1] += dd_buf[i][3] * coeff[0] + dd_buf[i][5] * coeff[1] * 0.5;
                                tosolve_den2.data[3] += dd_buf[i][4] * coeff[1] + dd_buf[i][5] * coeff[0] * 0.5;
                            }
                            d_scale[0][0] = tosolve_num2[1] + tosolve_den2.data[1] * curcol[mergesuggestion[0]];
                            d_scale[0][1] = tosolve_num2[2] + tosolve_den2.data[3] * curcol[mergesuggestion[0]];
                            candidates[0] = curcol;
                            candidates[0][mergesuggestion[1]] += curcol[mergesuggestion[0]];
                            candidates[0].erase(mergesuggestion[0]);
                        }
                    }
                    d_scale[0] = dd_scale * d_scale[0];
                    candidates[0].keysort();
                    sum[0] = 0.0f;
                    if (auto ite = candidates[0].getIterator()) {
                        do{
                            sum[0] += *ite;
                        }while(ite++);
                        if (fabs(sum[0] - 1.0) > 0.1) {printf("hidden state norm is hh %e!\n", sum[0]);
                            candidates[0].show();
                            curcol.show();
                        exit(1);}
                    }else {printf("got empty candidate!\n"); exit(1);}
                    double tmp;
                    //printf("cur coor %i/%i\n", coor[1], target.colhid.getNBcols());fflush(stdout);
                    lld[2] = 0.0f;
                    tosolve_num2[0] = 1.0 / curcol.getSize();
                    if (auto ite = curcol.getIterator()) do{
                        acoor[0]  = ite();
                        for(acoor[1]=0;acoor[1]< target.par_RM.sizes[1]; acoor[1]++){
                            coor[0] = acoor[1];
                            if (freq_colh[coor[0]] <= row_freq_hcol(coor)- 0.001) thrbase.msgs.insert(string("got negative!") + to_string(freq_colh[coor[0]] - row_freq_hcol(coor)));
                            lld[2] += tosolve_num2[0] * (row_freq_hcol(coor) * log(target.par_D(acoor)) + (freq_colh[coor[0]] - row_freq_hcol(coor)) * log(1.0 - target.par_D(acoor)));
                        }
                    }while(ite++);
                    lld[2] *= dropout_likelihood_factor;
                    if (!ExOp::isValid(lld[2])){
                        thrbase.msgs.insert(string("Default Dropout likelihood is NA"));
                    }

                    insist =0;
                    while(true){

                        // Step 4 Evaluate Quantidates
                        for(j=0;j<2;j++){
                            future_statistics[j].toZero();future_statistic_trix[j].toZero();
                            curproj[j] = target.par_RM.mkBackMult(candidates[j]);
                            lls[j] = target.row_states_prior(candidates[j]) * ((use_heavy_prior) ? trdata.data[coor[1]].getSize() : 1);
                            if (candidates[j].hasSameSparsity(curcol)) lld[j] = lld[2];
                            else{
                                lld[j] = 0.0f;
                                tosolve_num2[0] = 1.0 / candidates[j].getSize();
                                if (auto ite = candidates[j].getIterator()) do{
                                    acoor[0] = ite();
                                    for(acoor[1] =0;acoor[1]< target.par_RM.sizes[1]; acoor[1]++){
                                        coor[0] = acoor[1];

                                        if ((!ExOp::isValid(log(target.par_D(acoor))))||(!ExOp::isValid(log(1.0 - target.par_D(acoor))))) thrbase.msgs.insert(string("got problem! with " + to_string(acoor[0])+ string(" ")+  to_string(acoor[1]))) ;


                                        lld[j] += tosolve_num2[0] * (row_freq_hcol(coor) * log(target.par_D(acoor)) + (freq_colh[coor[0]] - row_freq_hcol(coor)) * log(1.0 - target.par_D(acoor)));
                                    }
                                }while(ite++);
                                lld[j] *= dropout_likelihood_factor;
                                if (!ExOp::isValid(lld[j])) thrbase.msgs.insert(string("Canditate dropoutLL is invalid! ")+ to_string(j));
                            }
                            if ((candidates[1].getSize() == 0)||(insist > 0)) break;
                        }
                                     //   printf("loopyb\n"); fflush(stdout);
                        for(i=0;i<trdata.data[coor[1]].getSize();i++){
                            coor[0]= trdata.data[coor[1]].deref_key(i);
                            if (!ExOp::isValid(target.par_colscale[coor[0]])) continue; // can do better...
                            for(j=0;j<2;j++){
                                // first check
                                coeff = target.par_rowscale[coor[1]];
                                coeff += d_scale[j];
                                if ((tmp = (coeff[0]*coeff[0] + coeff[1]*coeff[1])) > 100.0){ // too far, send to boundary
                                    tmp = pow(tmp, -0.5);
                                    d_scale[j][0] = coeff[0] * tmp;
                                    d_scale[j][1] = coeff[1] * tmp;
                                    d_scale[j] -= target.par_rowscale[coor[1]];
                                }

                                coeff = curproj[j].mkInnerProd(target.colhid.getColumn(coor[0]));
                                coeff += target.par_colscale[coor[0]];
                                coeff += target.par_rowscale[coor[1]];
                                coeff += d_scale[j];
                                deflatedNBterm_wr_f_d2(ld.evilterms, trdata.data[coor[1]].deref(i), coeff[0], coeff[1]);
                               if (auto ite = target.colhid.getColumn(coor[0]).getIterator()) do{
                                    future_statistics[j][ite()][0] += ld.evilterm_alt.d[0] * (*ite);
                                    future_statistics[j][ite()][1] += ld.evilterm_alt.d[1] * (*ite);

                                    trix_input[0] = ld.evilterms[3] * (*ite);
                                    trix_input[1] = ld.evilterms[4] * (*ite);
                                    trix_input[2] = ld.evilterms[5] * (*ite) * 0.5;
                                    if (auto ite2 = target.colhid.getColumn(coor[0]).getIterator()) do{
                                        if (ite2() > ite()) continue;
                                        future_statistic_trix[j].data[ite2() + ((ite() * (ite() + 1)) >> 1)] += trix_input * (*ite2);
                                    }while(ite2++);

                                }while(ite++);

                                lls[j] += ld.evilterms[0];
                                if ((candidates[1].getSize() == 0)||(insist > 0)) break;
                            }
                        }


                        // worrying about dropouts!
                        j = ((candidates[1].getSize() == 0)||(lls[0] + lld[0] > lls[1] + lld[1])) ? 0 : 1;
                        j = (lls[j] + lld[j] > lls[2] + lld[2]) ? j : 2;
                        if (j != 2){
                            if (candidates[j].hasSameSparsity(curcol)) update_type[coor[1]] = 'S';
                            else {
                                update_type[coor[1]]= 'D';
                                // registering change...
                                tosolve_num2[0] = 1.0 / curcol.getSize();
                                tosolve_num2[1] = 1.0 / candidates[j].getSize();
                                if (auto ite = curcol.getIterator()) do{ freq_diff[threadID][ite()] -= tosolve_num2[0];} while(ite++);
                                if (auto ite = candidates[j].getIterator()) do{acoor[0] = ite(); freq_diff[threadID][ite()] += tosolve_num2[1];} while(ite++);
                                for(i = 0; i < trdata.data[coor[1]].getSize();i++){
                                    tosolve_num2[0] = 1.0 / (curcol.getSize() *  target.colhid.data[trdata.data[coor[1]].deref_key(i)].getSize());
                                    tosolve_num2[1] = 1.0 / (candidates[j].getSize() *  target.colhid.data[trdata.data[coor[1]].deref_key(i)].getSize());
                                    if (auto ite2 = target.colhid.data[trdata.data[coor[1]].deref_key(i)].getIterator()) do{
                                        acoor[1] = ite2();
                                        if (auto ite = curcol.getIterator()) do{acoor[0] = ite(); nonzero_freq_diff[threadID](acoor) -= tosolve_num2[0];} while(ite++);
                                        if (auto ite = candidates[j].getIterator()) do{acoor[0] = ite(); nonzero_freq_diff[threadID](acoor) += tosolve_num2[1];} while(ite++);
                                    }while(ite2++);
                                }
                                nb_flags[threadID][2]++;
                            /*
                            nbmax_diff[threadID][candidates[j].deref_key(0)]++;
                            nbmax_diff[threadID][curcol.deref_key(0)]--;
                            for(i = 0; i < trdata.data[coor[1]].getSize();i++){
                                acoor[1] = target.colhid.data[trdata.data[coor[1]].deref_key(i)].deref_key(0);
                                acoor[0] = curcol.deref_key(0);
                                nonzero_diff[threadID](acoor)--;
                                acoor[0] = candidates[j].deref_key(0);
                                nonzero_diff[threadID](acoor)++;
                            }*/
                            }
                        }

                        if (j != 2) break;
                        if (insist == 0) {
                            if (candidates[1] == curcol) {thrbase.msgs.insert(string("Failed to update scale alone!\n")); update_type[coor[1]] = 'F'; break;}
                            candidates[0].toMemmove(candidates[1]);
                        }
                        if (insist++ > 10) {update_type[coor[1]] = 'F'; break;}

                        if (auto ite = curcol.getIterator()) do{
                            candidates[0][ite()] += *ite;
                        }while(ite++);

                        if (auto ite = candidates[0].getIterator()) do{
                            *ite *= 0.5;
                        }while(ite++);
                        d_scale[0] *= 0.5;
                    }
                    nb_flags[threadID][1] += insist;
            //printf("curss coor %i/%i\n", coor[1], target.colhid.getNBcols());fflush(stdout);

                    /*ll[threadID][0] -= target.row_states_prior(curcol);
                    ll[threadID][1] -= target.row_states_prior(candidates[j]);
                    ll[threadID][4] += target.row_states_prior(curcol);
                    ll[threadID][5] += target.row_states_prior(candidates[j]);*/

                    // Step 6 save scope for parameter update
                    newscale[coor[1]] = target.par_rowscale[coor[1]];
                    if (j == 2) { // failed to improve...
                        //printf("row %i failed %e %e\n", coor[1], target.par_rowscale[coor[1]][0], target.par_rowscale[coor[1]][1] );
                        candidates[0] = curcol;
                        newhidden.memmoveColumn(coor[1], candidates[0]);
                        nb_flags[threadID][0]++;
                    }else{
                        newhidden.memmoveColumn(coor[1], candidates[j]);
                        newscale[coor[1]] += d_scale[j];
                    }
                    ll[threadID][0] += lls[2];
                    ll[threadID][1] += lls[j];
                    ll[threadID][2] += lld[2];
                    ll[threadID][3] += lld[j];
                    ll[threadID][4] += lls[3];
                    candidates[1].toMemfree();

                    future_factor.toMemfree();
                    if (auto ite = newhidden.data[coor[1]].getIterator()) do{
                        future_factor[ite()][0] = *ite; future_factor[ite()][1] = *ite;
                    }while(ite++);
                    dll_dMR2[threadID].toAddBackOuterProd(future_statistics[j], future_factor);


                    uint32_t dashift = target.par_RM.sizes[0] * target.par_RM.sizes[1];
                    uint32_t k;
                    for(i = 0,k=0; i< target.par_RM.sizes[1];i++){
                        for(uint32_t l=0;l<=i;l++,k++){
                            if (auto ite = newhidden.data[coor[1]].getIterator()) do{
                                uint32_t m = ite() + i * target.par_RM.sizes[0];
                                double* currow_A = dd_massive[threadID].data + (((m * (m +1)) >> 1) + l * target.par_RM.sizes[0]);
                                double* currow_B = dd_massive[threadID].data + ((((m + dashift) * (m + dashift +1)) >> 1) + l * target.par_RM.sizes[0]);
                                if (auto ite2 = newhidden.data[coor[1]].getIterator()) do{
                                    if (ite2() > ite()) continue;
                                    double tmp = (*ite) * (*ite2);
                                    currow_A[ite2()] += future_statistic_trix[j].data[k][0] * tmp;
                                    currow_B[ite2()] += future_statistic_trix[j].data[k][2] * tmp;
                                    currow_B[ite2() + dashift] += future_statistic_trix[j].data[k][1] * tmp;
                                }while(ite2++);
                            }while(ite++);
                        }
                    }
                }
            }
            lls[0] =0;

            for(coor[1]=colrange[threadID];coor[1] < colrange[threadID+1];coor[1]++){
                if (data.getColSize(coor[1]) != 0) lls[0] += target.col_states_prior(target.colhid.getColumn(coor[1])) * ((use_heavy_prior) ? data.data[coor[1]].getSize() : 1);
			//log(target.colhid_prior.mkInnerMult(target.colhid.getColumn(coor[1]))+ target.colhid_prior_constant) * ((tmptmp) ? data.data[coor[1]].getSize() : 1);
            }
            ll[threadID][0] += lls[0];
            ll[threadID][1] += lls[0];
            ll[threadID][4] += lls[0];
        break;
        case 2:
            for(coor[1]=rowrange[threadID];coor[1] < rowrange[threadID+1];coor[1]++){ thrbase.updateProgress(threadID);
                if (auto ite = trdata.data[coor[1]].getIterator()) do{
                    if (update_type[ite()] == 'D'){
                        //msgs.insert(to_string(target.colhid.data[ite()].getSize()) + string("to") + to_string(newhidden.data[ite()].getSize()) );
                        tosolve_num2[0] = 1.0 / newhidden.data[ite()].getSize();
                        tosolve_num2[1] = 1.0 / target.colhid.data[ite()].getSize();
                        if (auto ite2 = newhidden.data[ite()].getIterator()) do {coor[0] = ite2(); row_freq_hcol(coor) += tosolve_num2[0];} while(ite2++);
                        if (auto ite2 = target.colhid.data[ite()].getIterator()) do {coor[0] = ite2(); row_freq_hcol(coor) -= tosolve_num2[1];} while(ite2++);
                    }
                }while(ite++);
            }
            /* to check the update...
            for(coor[1]=rowrange[threadID];coor[1] < rowrange[threadID+1];coor[1]++){ thrbase.updateProgress(threadID);
                if (auto ite = trdata.data[coor[1]].getIterator()) do{
                    tmp = 1.0 / newhidden.data[ite()].getSize();
                    if (auto ite2 = newhidden.data[ite()].getIterator()) do{coor[0] = ite2(); row_freq_hcol(coor) -= tmp;} while(ite2++);
                }while(ite++);
            }*/

        break;
        case 3:
            for(coor[1]=colrange[threadID];coor[1] < colrange[threadID+1];coor[1]++){ thrbase.updateProgress(threadID);
                if (auto ite = data.data[coor[1]].getIterator()) do{
                    if (update_type[ite()] == 'D'){
                    tosolve_num2[0] = 1.0 / newhidden.data[ite()].getSize();
                    tosolve_num2[1] = 1.0 / target.rowhid.data[ite()].getSize();
                    if (auto ite2 = newhidden.data[ite()].getIterator()) do {coor[0] = ite2(); col_freq_hrow(coor) += tosolve_num2[0];} while(ite2++);
                    if (auto ite2 = target.rowhid.data[ite()].getIterator()) do {coor[0] = ite2(); col_freq_hrow(coor) -= tosolve_num2[1];} while(ite2++);
                    }
                }while(ite++);
            }


        break;
        case 4:
            nonzero_freq_diff[threadID].toZero();
            for(i=colrange[threadID];i < colrange[threadID+1];i++){ thrbase.updateProgress(threadID);
                if (data.data[i].getSize() == 0){
                    if (target.colhid.data[i].getSize() != 0 ) {
                        printf("warning, had invalid state for empty collumn %i\n", i);
                        //target.colhid.data[i].toMemfree();
                    }
                }else{
                    if (auto ite = data.data[i].getIterator()) do{
                        coor[1] = i;
                        tmp = 1.0 / target.rowhid.data[ite()].getSize();
                        if (auto inrte = target.rowhid.data[ite()].getIterator()) do{coor[0] = inrte(); col_freq_hrow(coor) += tmp;} while(inrte++);
                        tmp = 1.0 / (target.colhid.data[i].getSize() * target.rowhid.data[ite()].getSize());
                        if (auto inrte = target.rowhid.data[ite()].getIterator()) do{
                            coor[0] = inrte();
                            if (auto inrte2 = target.colhid.data[i].getIterator()) do{coor[1] = inrte2(); nonzero_freq_diff[threadID](coor)+= tmp;} while(inrte2++);
                        }while(inrte++);
                    }while (ite++);
                }
            }

            for(coor[1]=rowrange[threadID];coor[1] < rowrange[threadID+1];coor[1]++){ thrbase.updateProgress(threadID);
                if (trdata.data[coor[1]].getSize() == 0){
                    if (target.rowhid.data[coor[1]].getSize() != 0 ) {
                        printf("warning, had invalid state for empty row %i\n", i);
                        //target.rowhid.data[coor[1]].toMemfree();
                    }
                }else{
                    if (auto ite = trdata.data[coor[1]].getIterator()) do{
                        tmp = 1.0 / target.colhid.data[ite()].getSize();
                        if (auto inrte = target.colhid.data[ite()].getIterator()) do{coor[0] = inrte(); row_freq_hcol(coor) += tmp;} while(inrte++);
                    }while(ite++);
                }
            }
        break;
        }
/*
        for(coor[1]=colrange[threadID];coor[1] < colrange[threadID+1];coor[1]++){
            acoor[1] = newhidden.data[coor[1]].deref_key(0);
            if (auto ite = data.data[coor[1]].getIterator()) do{
                acoor[0] = target.rowhid.data[ite()].deref_key(0);
                nonzero_diff_dbg[threadID](acoor)--;
            }while(ite++);
        }*/

        thrbase.finishProgress(threadID);
    }
};

LFHTEMP void BiClassifier<C,R,S>::run2D_EM_v4(const SparseMatrix<C>& data, const Vector<uint32_t> &excl_list,  uint32_t nbstep, bool use_heavy_prior){
    BiClassifier<C,R,S>::Task_run2D_EM_v4 lotask(*this,data,thrbase.nbthreads, use_heavy_prior);
    Tuple<uint32_t, 2u> coor, coor2, coorC, coorR;
    double tmp;
    uint32_t i,j,k;
    i=0;
    for(coor[1]=0;coor[1]<rowhid.getNBcols();coor[1]++) {
        if (lotask.trdata.data[coor[1]].getSize() == 0) continue;
        i += lotask.trdata.data[coor[1]].getSize();
        tmp = 1.0 / rowhid.data[coor[1]].getSize();
        if (auto ite = rowhid.data[coor[1]].getIterator()) do{lotask.freq_rowh[ite()] += tmp;}while(ite++);
    }
    // the matrix has rowhid.getNBcols() x colhid.getNBcols() observations, but only "i" data points, since we do not
    // want to focus too much on dropouts, we define "dropout_likelihood_factor" as weight for the dropout likelihoods
    lotask.dropout_likelihood_factor = (((double)i) / rowhid.getNBcols() ) / colhid.getNBcols();
    printf("dropout likelihood weight: %e\n", lotask.dropout_likelihood_factor);
    for(coor[1]=0;coor[1]<colhid.getNBcols();coor[1]++){
        if (lotask.data.data[coor[1]].getSize() == 0) continue;
        tmp = 1.0 / colhid.data[coor[1]].getSize();
        if (auto ite = colhid.data[coor[1]].getIterator()) do{lotask.freq_colh[ite()] += tmp;}while(ite++);
    }




    lotask.taskID = 4;
    thrbase.startProgress("Populating Search Scope", lotask.colrange.last() + lotask.rowrange.last());
    for(i=thrbase.nbthreads-1;i>0;i--) thrbase.submit(lotask , i);
    thrbase.submit_ThenWait(lotask, 0);
/*
    printf("that's the start...row\n");
    for(coor[1]=0;coor[1] < 40; coor[1]++) {
        tmp =0;
        for(coor[0]=0;coor[0] < par_RM.sizes[1]; coor[0]++) tmp += lotask.row_freq_hcol(coor);
        printf("%i:%e (size is %i)\n", coor[1], tmp,lotask.trdata.data[coor[1]].getSize());
    }
    printf("end..colw\n");
    for(coor[1]=0;coor[1] < 40; coor[1]++) {
        tmp =0;
        for(coor[0]=0;coor[0] < par_RM.sizes[0]; coor[0]++) tmp += lotask.col_freq_hrow(coor);
        printf("%i:%e (size is %i)\n", coor[1], tmp,lotask.data.data[coor[1]].getSize());
    }*/

    /*
    for(coor[1]=0;coor[1] < 200;coor[1]++) printf("%i:%e\n", lotask.row_freq_hcol.data[coor[1]]);
    if (auto ite = lotask.data.data[0].getIterator()) do{printf("C0, %i:%i\n",ite(),*ite);}while(ite++);
    if (auto ite = lotask.data.data[1].getIterator()) do{printf("C1, %i:%i\n",ite(),*ite);}while(ite++);
    if (auto ite = lotask.data.data[2].getIterator()) do{printf("C2, %i:%i\n",ite(),*ite);}while(ite++);
    if (auto ite = lotask.trdata.data[0].getIterator()) do{printf("R0, %i:%i\n",ite(),*ite);}while(ite++);
    if (auto ite = lotask.trdata.data[1].getIterator()) do{printf("R1, %i:%i\n",ite(),*ite);}while(ite++);
    if (auto ite = lotask.trdata.data[2].getIterator()) do{printf("R2, %i:%i\n",ite(),*ite);}while(ite++);*/
   /* if (auto ite = colhid.data[0].getIterator()) do{printf("H0, %i:%e\n",ite(),*ite);}while(ite++);
    if (auto ite = colhid.data[1].getIterator()) do{printf("H1, %i:%e\n",ite(),*ite);}while(ite++);
    if (auto ite = colhid.data[2].getIterator()) do{printf("H2, %i:%e\n",ite(),*ite);}while(ite++);
    if (auto ite = rowhid.data[0].getIterator()) do{printf("Z0, %i:%e\n",ite(),*ite);}while(ite++);
    if (auto ite = rowhid.data[1].getIterator()) do{printf("Z1, %i:%e\n",ite(),*ite);}while(ite++);
    if (auto ite = rowhid.data[2].getIterator()) do{printf("Z2, %i:%e\n",ite(),*ite);}while(ite++);*/

    lotask.nonzero_freq_grid = lotask.nonzero_freq_diff[0];
    for(i=1;i<lotask.alphas.getSize();i++) {
        lotask.nonzero_freq_grid += lotask.nonzero_freq_diff[i];
    }
/*
     for(int i=0;i < data.data.getSize();i++){
        coor[1] = colhid.data[i].deref_key(0);
        if (auto ite = data.data[i].getIterator()) do{
            if (colhid.data[i].getSize() == 0) printf("impossible... col?\n");
            if (rowhid.data[ite()].getSize() == 0) printf("impossible... row?\n");
            coor[0] = colhid.data[i].deref_key(0);
            coor[1] = ite();
            //lotask.row_nbmax_hcol(coor)++;

            coor[0] = rowhid.data[ite()].deref_key(0);
            coor[1] = i;
          //  lotask.col_nbmax_hrow(coor)++;


            coor[1] = colhid.data[i].deref_key(0);
            lotask.maxgrid_nonzero_count(coor)++;
        }while (ite++);
    }*/
    /*tmp = 0; if (auto ite = lotask.freq_rowh.getIterator()) do{tmp += *ite;} while(ite++);
    printf("freq_rowh sum:%f\n", tmp);
    lotask.freq_rowh.show();
    tmp = 0; if (auto ite = lotask.freq_colh.getIterator()) do{tmp += *ite;} while(ite++);
    printf("freq_colh sum:%f\n", tmp);
    lotask.freq_colh.show();
    tmp = 0; if (auto ite = lotask.nonzero_freq_grid.getIterator()) do{tmp += *ite;} while(ite++); printf("nonzero_freq_grid sum:%f\n", tmp);
    */

    double last_dll = 0.0;
    double dropLL;


    if (auto ite = par_D.getIterator()) do{
        tmp = ((double)(lotask.freq_rowh[ite()[0]] * lotask.freq_colh[ite()[1]]));
        (*ite) = (lotask.nonzero_freq_grid(ite())+0.5) / (tmp + 1.0);
        double  term = log(lotask.nonzero_freq_grid(ite())) * lotask.nonzero_freq_grid(ite());
                term += log(tmp - lotask.nonzero_freq_grid(ite())) * (tmp - lotask.nonzero_freq_grid(ite()));
                term -= log(tmp) * tmp;
        //printf("%e %e %e is the rest...\n", lotask.nonzero_freq_grid(ite()), tmp - lotask.nonzero_freq_grid(ite()), tmp);
        if (ExOp::isValid(term)) last_dll += term; // value is invalid if there is no entropy for that position, LL is 0 if so
    }while(ite++);
    printf("Original dropLL %e\n", last_dll * lotask.dropout_likelihood_factor);
    double that_ll_constant=0.0;
    if (auto ite = data.getIterator()) do{
        that_ll_constant -= lngamma(1.0 + *ite);
    }while(ite++);


    lotask.ignore_drop = false;

    coor[0] = par_RM.sizes[0];
    coor[1] = par_RM.sizes[1];

    double old_LL = -DBL_MAX;


    double tmpsum;
    double ll,ll_best;
    FunctionScope<double&,double&,double> funabsoft(ExCo<double>::toAbsoftInverse,pow(0.5, 256.0f));
    Trianglix<double,2u> single_denum;
    Tuple<double, 2u> single_num;
    Tuple<double,0u,TUPLE_FLAG_REMOTE_MEMORY> dainput;
    uint32_t tmpsize;
    int step;
    CurvatureSearchScope cs;
    uint32_t fuz =0;
    double expected =0.0;
    double expected2 =0.0;
    Tuple<uint32_t> permute;

    double rectsize;

    Tuple<double,0u> toproj; toproj.setSize(2 * par_RM.sizes[0] * par_RM.sizes[1]);
    cs.init(2 * par_RM.sizes[0] * par_RM.sizes[1], 1.0);


    Tuple<uint32_t> spcount;

    double logalpha = -1.0;

    lotask.taskID =0;
    double dadada;
    lotask.oldRM = par_RM;
    while(true){

        #ifdef Rcpp_hpp
            Rcpp::checkUserInterrupt();
        #endif // Rcpp_hpp


        lotask.newhidden.setNBcols((lotask.taskID == 0)? colhid.getNBcols() : rowhid.getNBcols());
        lotask.newscale.setSize(lotask.newhidden.getNBcols());
 //       scp.new_nbmax_h.setSizes((scp.istimefor_columns)? scp.row_nbmax_hcol.getDims() : scp.col_nbmax_hrow.getDims()  );

        lotask.tmptmp = false; //(fuz >= (nbstep/2)); // heavy/light hidden prior
        //scp.heuristictime = true;

        if ((fuz % 10) == 5) doRecenter(); // hehe

       /* if (fuz != 0){
        lotask.taskID ^= 1;
        thrbase.startProgress(((lotask.taskID == 0) ? "Fake Updating Columns" : "Fake Updating Rows"), (lotask.taskID == 0) ? lotask.colrange.last() : lotask.rowrange.last() );

        //for(i=thrbase.nbthreads-1;i>0;i--) thrbase.submit(std::bind(&Task::operator(), std::ref(lotask), i) );
        //thrbase.submit_ThenWait(std::bind(&Task::operator(), std::ref(lotask), 0));

        //thrbase.submit_Array(lotask,thrbase.nbthreads);
        for(i=thrbase.nbthreads-1;i>0;i--) thrbase.submit(lotask , i);
        thrbase.submit_ThenWait(lotask, 0);
        for(i=1;i<lotask.alphas.getSize();i++) lotask.ll[0] += lotask.ll[i];
        printf("LL is %e  (base %e), (d-out) %e,   %e old old\n",lotask.ll[0][2] + lotask.ll[0][0]  +that_ll_constant, lotask.ll[0][0] + that_ll_constant, lotask.ll[0][2], lotask.ll[0][4] + that_ll_constant);
        lotask.taskID ^= 1;
        }*/

        thrbase.startProgress(((lotask.taskID == 0) ? "Updating Columns" : "Updating Rows"), (lotask.taskID == 0) ? lotask.colrange.last() : lotask.rowrange.last() );

        //for(i=thrbase.nbthreads-1;i>0;i--) thrbase.submit(std::bind(&Task::operator(), std::ref(lotask), i) );
        //thrbase.submit_ThenWait(std::bind(&Task::operator(), std::ref(lotask), 0));

        //thrbase.submit_Array(lotask,thrbase.nbthreads);
        for(i=thrbase.nbthreads-1;i>0;i--) thrbase.submit(lotask , i);
        thrbase.submit_ThenWait(lotask, 0);

       // (lotask.nonzero_diff[0] + lotask.nonzero_diff_dbg[0]).show();
        tmpsize = (lotask.taskID == 0) ? par_RM.sizes[1]: par_RM.sizes[0];
        for(i=1;i<lotask.alphas.getSize();i++) {
            lotask.ll[0] += lotask.ll[i];
            lotask.alphas[0] += lotask.alphas[i];
            lotask.nonzero_freq_diff[0] += lotask.nonzero_freq_diff[i];
            lotask.freq_diff[0] += lotask.freq_diff[i];
            lotask.nb_flags[0] += lotask.nb_flags[i];
        }
      //  printf("total local ll %e %e\n", lotask.ll[0][0], lotask.ll[0][1]);
        if (!lotask.ignore_drop) {
            //printf("cached change\n");
            //lotask.nonzero_freq_diff[0].show();
            lotask.nonzero_freq_diff[0] += lotask.nonzero_freq_grid;

            dropLL =0.0;
            if (lotask.taskID == 0){
                for(i=0;i<par_RM.sizes[1];i++) lotask.freq_diff[0][i] += lotask.freq_colh[i];
                if (auto ite = par_D.getIterator()) do{
                    tmp = ((lotask.freq_rowh[ite()[0]] * lotask.freq_diff[0][ite()[1]]));
                    double term = log(*ite) * lotask.nonzero_freq_diff[0](ite());
                            term += log(1.0 - *ite) * (tmp - lotask.nonzero_freq_diff[0](ite()));
                    if (ExOp::isValid(term)) dropLL += term; // value is invalid if there is no entropy for that position, LL is 0 if so
                }while(ite++);
            }else{
                for(i=0;i<par_RM.sizes[0];i++) lotask.freq_diff[0][i] += lotask.freq_rowh[i];
                if (auto ite = par_D.getIterator()) do{
                    tmp = ((double)(lotask.freq_diff[0][ite()[0]] * lotask.freq_colh[ite()[1]]));
                    double  term = log(*ite) * lotask.nonzero_freq_diff[0](ite());
                            term += log(1.0 - *ite) * (tmp - lotask.nonzero_freq_diff[0](ite()));
                    if (ExOp::isValid(term)) dropLL += term; // value is invalid if there is no entropy for that position, LL is 0 if so
                }while(ite++);
            }
            last_dll = 0.0;
            if (auto ite = par_D.getIterator()) do{
                tmp = ((lotask.freq_rowh[ite()[0]] * lotask.freq_colh[ite()[1]]));
                double term = log(*ite) * lotask.nonzero_freq_grid(ite());
                        term += log(1.0 - *ite) * (tmp - lotask.nonzero_freq_grid(ite()));
                if (ExOp::isValid(term)) last_dll += term; // value is invalid if there is no entropy for that position, LL is 0 if so
            }while(ite++);
            last_dll *= lotask.dropout_likelihood_factor;
            dropLL *= lotask.dropout_likelihood_factor;

            printf("the LL %e -> %e;\n", lotask.ll[0][0]+ that_ll_constant, lotask.ll[0][1]+ that_ll_constant);
            printf("the dropout ll %e -> %e\n", lotask.ll[0][2], lotask.ll[0][3]);
            printf("chk dropout ll %e -> %e\n", last_dll, dropLL);
            dadada = lotask.ll[0][1] + that_ll_constant;
            lotask.ll[0][0] += lotask.ll[0][2]; // not exact (dropout), but not important
            lotask.ll[0][1] += lotask.ll[0][3];
        }
        lotask.ll[0][0] += that_ll_constant;
        lotask.ll[0][1] += that_ll_constant;
        lotask.ll[0][4] += that_ll_constant;
        printf("total local ll %e %e\n", lotask.ll[0][0], lotask.ll[0][1]);

        lotask.alphas[0] /= lotask.alphas.getSize();
        for(i=1;i<lotask.alphas.getSize();i++) lotask.alphas[i] = lotask.alphas[0] ;

        rectsize = ((double)lotask.nb_flags[0][0]) / ((lotask.taskID == 0) ? lotask.data.getNBcols() : lotask.trdata.getNBcols());
        rectsize = ((lotask.taskID == 0) ? lotask.data.getNBcols() : lotask.trdata.getNBcols());
        printf("Fraction of Sparse Move %e Local move %e No move %e\n", ((double)lotask.nb_flags[0][2]) / rectsize, 1.0 - ((double)lotask.nb_flags[0][0]+lotask.nb_flags[0][2]) / rectsize, ((double)lotask.nb_flags[0][0]) / rectsize);
        printf("da old LL should be %e (excludes dropout LL)\n", lotask.ll[0][4]);
        printf("%i Likelihood %e -> %e -> %e (%c In rate, insist %e)\n", fuz, old_LL, lotask.ll[0][0], lotask.ll[0][1], (lotask.taskID == 0) ? 'C' : 'R', ((double)lotask.nb_flags[0][1]) / rectsize);
        printf("expected change %e, got %e  (ddmag %e)\n", expected + expected2, lotask.ll[0][0] - old_LL, expected2);

        if ((old_LL < lotask.ll[0][1])||(fuz == 0)||((lotask.tmptmp)&&(fuz == (nbstep/2)))){
            expected = 0.0;
            expected2 = 0.0;
            lotask.oldRM = par_RM;
            old_LL = lotask.ll[0][1];

            lotask.taskID |= 2;


            thrbase.startProgress(((lotask.taskID == 2) ? "Accepting new column states" : "Accepting new row states"), (lotask.taskID == 2) ? lotask.rowrange.last() : lotask.colrange.last());
            for(i=thrbase.nbthreads-1;i>0;i--) thrbase.submit(lotask , i);
            thrbase.submit_ThenWait(lotask, 0);

          //  lotask.freq_diff[0].show();
          //  printf("dachk\n");
          //  for(coor[1]=0;coor[1] < 20; coor[1]++) for(coor[0]=0;coor[0] < par_RM.sizes[1]; coor[0]++) printf("%i,%i:%e (d=%e)\n", coor[0], coor[1], lotask.row_freq_hcol(coor), lotask.freq_diff[0][coor[0]]- lotask.row_freq_hcol(coor));

             /*/lotask.taskID = 4;
             colhid.toMemmove(lotask.newhidden);
            thrbase.startProgress(((lotask.taskID == 2) ? "Accepting new column states" : "Accepting new row states"), (lotask.taskID == 2) ? lotask.rowrange.last() : lotask.colrange.last());
            for(i=thrbase.nbthreads-1;i>0;i--) thrbase.submit(lotask , i);
            thrbase.submit_ThenWait(lotask, 0);*/

            //lotask.taskID = 2;
            if (lotask.taskID == 2){
/*
                    printf("and dachk\n");
                    for(coor[1]=1700;coor[1] < 1780; coor[1]++) {
                        tmp =0;
                        for(coor[0]=0;coor[0] < par_RM.sizes[1]; coor[0]++) tmp += lotask.row_freq_hcol(coor);
                        printf("%i:%e (size is %i)\n", coor[1], tmp,lotask.trdata.data[coor[1]].getSize());
                    }*/

                    /*
                for(i=0;i< colhid.getNBcols();i++){
                    if ((colhid.data[i].getSize() == 0)||(lotask.newhidden.data[i].getSize() == 0)) continue;
                    if (colhid.data[i].deref_key(0) != lotask.newhidden.data[i].deref_key(0)){
                        for(j = 0; j < data.data[i].getSize();j++){
                            coor[1] = data.data[i].deref_key(j);
                            coor[0] = colhid.data[i].deref_key(0);
                            lotask.row_nbmax_hcol(coor)--;
                            coor[0] = lotask.newhidden.data[i].deref_key(0);
                            lotask.row_nbmax_hcol(coor)++;
                        }
                    }
                }*/

                /*uint32_t nbnb =0;
                if (auto ite = lotask.row_freq_hcol.getIterator()) do{
                    if ((!ExOp::isValid(*ite))||(*ite) < 0) nbnb++;
                }while(ite++);
                printf("number of problems %i\n",nbnb);*/


                for(i=0;i<par_RM.sizes[1];i++) lotask.freq_colh[i] = lotask.freq_diff[0][i];
                colhid.toMemmove(lotask.newhidden);
                par_colscale.toMemmove(lotask.newscale);
                permute = col_states_prior.learnPermutation(colhid,1.0, dadada);
            //    permute.setSize(10);
            //    for(int iii=0;iii<10;iii++) permute[iii] = iii;




                /*
                auto lscp = col_states_prior.mkEMscope();
                lscp.init();
                for(coor[1]=0;coor[1] < colhid.getNBcols();coor[1]++) {
                    if (colhid.data[coor[1]].getSize() != 0) lscp.EMregist(colhid.getColumn(coor[1]));
                }
                col_states_prior.learn(lscp,1.0);*/


             }else{
                    /*
                    printf("end..colw\n");
                    for(coor[1]=1000;coor[1] < 1080; coor[1]++) {
                        tmp =0;
                        for(coor[0]=0;coor[0] < par_RM.sizes[0]; coor[0]++) tmp += lotask.col_freq_hrow(coor);
                        printf("%i:%e (size is %i)\n", coor[1], tmp,lotask.data.data[coor[1]].getSize());
                    }*/

                 /*
                for(i=0;i< rowhid.getNBcols();i++){
                    if ((rowhid.data[i].getSize() == 0)||(lotask.newhidden.data[i].getSize() == 0)) continue;
                    if (rowhid.data[i].deref_key(0) != lotask.newhidden.data[i].deref_key(0)){
                        for(j = 0; j < lotask.trdata.data[i].getSize();j++){
                            coor[1] = lotask.trdata.data[i].deref_key(j);
                            coor[0] = rowhid.data[i].deref_key(0);
                            lotask.col_nbmax_hrow(coor)--;
                            coor[0] = lotask.newhidden.data[i].deref_key(0);
                            lotask.col_nbmax_hrow(coor)++;
                        }
                    }
                }*/
                for(i=0;i<par_RM.sizes[0];i++) lotask.freq_rowh[i] = lotask.freq_diff[0][i];
                rowhid.toMemmove(lotask.newhidden);
                par_rowscale.toMemmove(lotask.newscale);
                permute = row_states_prior.learnPermutation(rowhid,1.0, dadada);

                /*auto lscp = row_states_prior.mkEMscope();
                lscp.init();
                for(coor[1]=0;coor[1] < rowhid.getNBcols();coor[1]++) {
                    if (rowhid.data[coor[1]].getSize() != 0) lscp.EMregist(rowhid.getColumn(coor[1]));
                }
                row_states_prior.learn(lscp,1.0); //printf("new col state likelihood %e\n", ); // be warned inside
                */

            }



            for(i=1;i<lotask.alphas.getSize();i++) {
                lotask.dd_massive[0] += lotask.dd_massive[i];
                lotask.dll_dMR2[0] += lotask.dll_dMR2[i];
            }
            // uptate D, close close form hurray!


            coor[0] = par_RM.sizes[0];
            coor[1] = par_RM.sizes[1];

           /* TMatrix<double> lochk;
            lochk.setSizes(coor).toZero();

            for(i=0;i <colhid.getNBcols();i++){
                if (auto ite = data.data[i].getIterator()) do{
                    coor[1] = i;
                    tmp = 1.0 / (colhid.data[i].getSize() * rowhid.data[ite()].getSize());
                    if (auto inrte = rowhid.data[ite()].getIterator()) do{
                        coor[0] = inrte();
                        if (auto inrte2 = colhid.data[i].getIterator()) do{coor[1] = inrte2(); lochk(coor) += tmp;} while(inrte2++);
                    }while(inrte++);
                }while (ite++);
            }
            printf("actual change\n");*/
//            (lochk - lotask.nonzero_freq_grid).show();

            lotask.nonzero_freq_grid = lotask.nonzero_freq_diff[0];
//            lotask.nonzero_freq_grid = lochk;

           /* tmp = 0; if (auto ite = lotask.freq_rowh.getIterator()) do{tmp += *ite;} while(ite++);
            printf("freq_rowh sum:%f\n", tmp);
            lotask.freq_rowh.show();
            tmp = 0; if (auto ite = lotask.freq_colh.getIterator()) do{tmp += *ite;} while(ite++);
            printf("freq_colh sum:%f\n", tmp);
            lotask.freq_colh.show();
            tmp = 0; if (auto ite = lotask.nonzero_freq_grid.getIterator()) do{tmp += *ite;} while(ite++); printf("nonzero_freq_grid sum:%f\n", tmp);

            lotask.freq_rowh.show();
            lotask.freq_colh.show();
            lotask.freq_rowh.toZero();*/
            /*
            for(coor[1]=0;coor[1]<rowhid.getNBcols();coor[1]++) {
                tmp = 1.0 / rowhid.data[coor[1]].getSize();
                if (auto ite = rowhid.data[coor[1]].getIterator()) do{lotask.freq_rowh[ite()] += tmp;}while(ite++);

            }
            lotask.freq_colh.toZero();
            for(coor[1]=0;coor[1]<colhid.getNBcols();coor[1]++){
                tmp = 1.0 / colhid.data[coor[1]].getSize();
                if (auto ite = colhid.data[coor[1]].getIterator()) do{lotask.freq_colh[ite()] += tmp;}while(ite++);
            }*/
            //printf("Chk chke!!!\n");
            //lotask.freq_rowh.show();
            //lotask.freq_colh.show();
            //lotask.nonzero_freq_grid.show();


            last_dll = 0.0;
            if (auto ite = par_D.getIterator()) do{
                (*ite) = ((double)lotask.nonzero_freq_grid(ite()) + 0.5) / (lotask.freq_rowh[ite()[0]] * lotask.freq_colh[ite()[1]] + 1.0);
                tmp = ((lotask.freq_rowh[ite()[0]] * lotask.freq_colh[ite()[1]]));
                double term = log(*ite) * lotask.nonzero_freq_grid(ite());
                        term += log(1.0 - *ite) * (tmp - lotask.nonzero_freq_grid(ite()));
                if (ExOp::isValid(term)) last_dll += term; // value is invalid if there is no entropy for that position, LL is 0 if so
            }while(ite++);
            printf("the revised ==>  %e  and %e  <==  bringing to %e\n", dadada, lotask.dropout_likelihood_factor * last_dll , dadada + lotask.dropout_likelihood_factor * last_dll);


            //printf("new dropout matrix:\n");
            //par_D.show();

            if (fuz >= nbstep) break;

            // fill symmetries... yea
            if (auto ite = lotask.dd_massive[0].getIterator()) do{
                if ((ite()[0] % par_RM.sizes[0]) > (ite()[1] % par_RM.sizes[0])){
                    (*ite) = lotask.dd_massive[0].cell( ((ite()[1] % par_RM.sizes[0]) + ((ite()[0] / par_RM.sizes[0]) * par_RM.sizes[0])) , ((ite()[0] % par_RM.sizes[0]) + ((ite()[1] / par_RM.sizes[0]) * par_RM.sizes[0])) );
                }
            }while(ite++);
            if (auto ite = lotask.dd_massive[0].getIterator()) do{
                if ((ite()[0] % (par_RM.sizes[0] * par_RM.sizes[1])) > (ite()[1] % (par_RM.sizes[0] * par_RM.sizes[1]))){
                    (*ite) = lotask.dd_massive[0].cell(ite()[1] % (par_RM.sizes[0] * par_RM.sizes[1]), (par_RM.sizes[0] * par_RM.sizes[1]) + ite()[0] % (par_RM.sizes[0] * par_RM.sizes[1]));
                }
            }while(ite++);

            /*
            i=0;
            if (auto ite = parascope[0].dll_dMR2.getIterator()) do{
                toproj[i++] = (*ite)[0];
                toproj[i++] = (*ite)[1];
            }while(ite++);
            cs.updateAscent(parascope[0].ll[1] - dropLL, &(par_RM.data[0][0]), &(toproj[0]));
            */

            i=0;
            if (auto ite = lotask.dll_dMR2[0].getIterator()) do{
                toproj[i + par_RM.sizes[0] * par_RM.sizes[1]] = (*ite)[1];
                toproj[i++] = (*ite)[0];
            }while(ite++);
            lotask.dd_massive[0].toEigenTransform(funabsoft);
            toproj = lotask.dd_massive[0] * toproj;


            #ifdef Rcpp_hpp
          //  if (lotask.taskID == 2){
            i=0;
                if (auto ite = par_RM.getIterator()) do{
                    expected += exp(logalpha) * toproj[i] * lotask.dll_dMR2[0](ite())[0];
                    expected += exp(logalpha) * toproj[i + par_RM.sizes[0] * par_RM.sizes[1]] * lotask.dll_dMR2[0](ite())[1];
                    (*ite)[1] += exp(logalpha) * toproj[i + par_RM.sizes[0] * par_RM.sizes[1]];
                    (*ite)[0] += exp(logalpha) * toproj[i++];
                }while(ite++);
           // }

            #endif

            /*

            if (auto ite = par_RM.getIterator()) do{
                single_num[0] = parascope[0].dll_dMR2(ite())[0];
                single_num[1] = parascope[0].dll_dMR2(ite())[1];
                single_denum.data[0] = parascope[0].dll_dMR2(ite())[2];
                single_denum.data[1] = parascope[0].dll_dMR2(ite())[4];
                single_denum.data[2] = parascope[0].dll_dMR2(ite())[3];
                single_denum.toEigenTransform(funabsoft);
                single_num = single_denum * single_num;

                if (ExOp::isValid(single_num)){
                    single_num *= 0.1;
                    //single_num[1] = (single_num[1] + single_num[0]);
                    //single_num[0] = single_num[1];
                    if ((*ite)[0] + single_num[0] < -10.0f) single_num[0] = -10.0 - (*ite)[0]; // dont let R too far...

                    (*ite) += single_num;
                    expected += single_num[0] * parascope[0].dll_dMR2(ite())[0];
                    expected += single_num[1] * parascope[0].dll_dMR2(ite())[1];
                    expected2 += single_num[0] * single_num[0] * parascope[0].dll_dMR2(ite())[2];
                    expected2 += single_num[1] * single_num[1] * parascope[0].dll_dMR2(ite())[3];
                    expected2 += 2.0 * single_num[0] * single_num[1] * parascope[0].dll_dMR2(ite())[4];

                }else{
                   printf("produced nans from:\n");
                   parascope[0].dll_dMR2(ite()).show();
                }

            }while(ite++);*/
            //expected += expected2;
            if (lotask.taskID == 2){
                colhid.permuteRows(permute);
                par_RM.permuteCols(permute);
                lotask.oldRM.permuteCols(permute);

                par_D.permuteCols(permute);
                lotask.nonzero_freq_grid.permuteCols(permute);
                lotask.row_freq_hcol.permuteRows(permute);
                lotask.freq_colh.permute(permute);
            }else{
                rowhid.permuteRows(permute);
                par_RM.permuteRows(permute);
                lotask.oldRM.permuteRows(permute);
                par_D.permuteRows(permute);
                lotask.nonzero_freq_grid.permuteRows(permute);
                lotask.col_freq_hrow.permuteRows(permute);
                lotask.freq_rowh.permute(permute);
            }


            lotask.taskID = 3^ lotask.taskID;
            logalpha = (logalpha * 127.0) / 128.0;

        }else{ // failed... fall backward and do not update
            if (fuz >= nbstep + 10) break;
            printf("was tried:\n");
            par_RM.show();
            printf("backup: %e\n", logalpha);
            lotask.oldRM.show();
            printf("FAILFAILFAILFAILFAILFAILFAILFAILFAILFAILFAIL %e\n", old_LL + that_ll_constant); fflush(stdout);
            par_RM = (lotask.oldRM * 3.0 + par_RM) * 0.25;
            printf("FAILFAILFAILFAILFAILFAILFAILFAILFAILFAILFAIL %e\n", lotask.ll[0][0] + that_ll_constant); fflush(stdout);
            expected *= 0.25;
            logalpha -= 1.0;
        }
        fuz++;
     }



     if (lotask.ignore_drop){ // update this, was not used
         /*lotask.maxgrid_nonzero_count.toZero();
         if (auto ite = data.getIterator()) do{
            coor[0] = rowhid.data[ite()[0]].deref_key(0);
            coor[1] = colhid.data[ite()[1]].deref_key(0);
            lotask.maxgrid_nonzero_count(coor)++;
         }while (ite++);
         lotask.nb_maxrowh.toZero();
         lotask.nb_maxcolh.toZero();
         for(coor[1]=0;coor[1]<rowhid.getNBcols();coor[1]++) {
            rowhid.makeMaximumAsFirst(coor[1]);
            lotask.nb_maxrowh[rowhid.data[coor[1]].deref_key(0)]++;
        }
        for(coor[1]=0;coor[1]<colhid.getNBcols();coor[1]++){
            colhid.makeMaximumAsFirst(coor[1]);
            lotask.nb_maxcolh[colhid.data[coor[1]].deref_key(0)]++;
        }
         if (auto ite = par_D.getIterator()) do{
            (*ite) = ((double)lotask.maxgrid_nonzero_count(ite())) / (lotask.nb_maxrowh[ite()[0]] * lotask.nb_maxcolh[ite()[1]]);
         }while(ite++);*/
    }




}


LFHTEMP void BiClassifier<C,R,S>::run2D_EM_v3(ThreadBase &tb, const SparseMatrix<C>& data, const Vector<uint32_t> &excl_list, const TMatrix<double>& soup, const Vector<uint32_t> &soupindex,  uint32_t nbstep, unsigned int nbthread){ class Task{
    public:
    BiClassifier<C,R,S>& target;
    TMatrix<uint32_t> maxgrid_nonzero_count; // (par_RM.sizes[0])
    DataGrid<uint32_t> row_nbmax_hcol; // (par_RM.sizes[1] x data.size[0])
    DataGrid<uint32_t> col_nbmax_hrow; // (par_RM.sizes[0] x data.size[1])
    Tuple<uint32_t> nb_maxrowh; //(par_RM.sizes[0])
    Tuple<uint32_t> nb_maxcolh; //(par_RM.sizes[1])

    Tuple< Tuple<double ,2u> > center_shift;

    SparseMatrix<double> introns;
    SparseMatrix<double> data;
    SparseMatrix<double> trdata;

    Matrix<double> bg_data;
    Tuple<uint32_t> bg_index;


    bool istimefor_columns;
    bool tmptmp;

    bool ignore_drop;
    double hiddenalpha;
    //bool heuristictime;

    Tuple<uint32_t> colrange;
    Tuple<uint32_t> rowrange;
    Tuple<Parallel_EMSCOPE> parascope;

    SparseMatrix<double> newhidden; // shared output
    Tuple< Tuple<double ,2u> > newscale; // shared output
    Task(BiClassifier<C,R,S>& _target, const SparseMatrix<C>& _data,const Vector<uint32_t> &excl_list):target(_target){
        if (excl_list.getSize() == 0){
        data = _data;
    }else{
        data.setNBcols(_data.getNBcols());
        myHashmap<uint32_t, void> which;
        for(uint32_t i=0;i<excl_list.getSize();i++) which.addEntry(excl_list[i]);
        for(uint32_t i=0;i<_data.getNBcols();i++){
            if (auto ite = _data.data[i].getIterator()) do{
                if (which.find(ite()) == 0xFFFFFFFF) data.data[i][ite()] = *ite;
            } while(ite++);
        }
    }
    trdata = data.mkTranspose();
    printf("fixing %i and %i nbrows so they match\n", trdata.data.getSize(), target.rowhid.getNBcols());
    while(trdata.data.getSize() < target.rowhid.getNBcols()) trdata.data.push_back();


    Tuple<uint32_t, 2u> coor;
    coor[0] = _target.par_RM.sizes[0];
    coor[1] = data.getNBcols();
    col_nbmax_hrow.toSizes(coor).toZero();

    coor[0] = _target.par_RM.sizes[1];
    coor[1] = trdata.getNBcols();
    row_nbmax_hcol.toSizes(coor).toZero();

    coor[0] = _target.par_RM.sizes[0];
    coor[1] = _target.par_RM.sizes[1];
    maxgrid_nonzero_count.setSizes(coor).toZero();

    nb_maxrowh.setSize(_target.par_RM.sizes[0]).toZero();
    nb_maxcolh.setSize(_target.par_RM.sizes[1]).toZero();

    center_shift.setSize( (_target.par_RM.sizes[0]) > _target.par_RM.sizes[1] ? _target.par_RM.sizes[0] : _target.par_RM.sizes[1] ).toZero();

    istimefor_columns = true;

    uint32_t i,j;
    for(coor[1]=0;coor[1]<_target.rowhid.getNBcols();coor[1]++) {
        _target.rowhid.makeMaximumAsFirst(coor[1]);
        nb_maxrowh[_target.rowhid.data[coor[1]].deref_key(0)]++;
    }
    for(coor[1]=0;coor[1]<_target.colhid.getNBcols();coor[1]++){
        _target.colhid.makeMaximumAsFirst(coor[1]);
        nb_maxcolh[_target.colhid.data[coor[1]].deref_key(0)]++;
    }

    if (auto ite = data.getIterator()) do{

        coor[0] = _target.colhid.data[ite.getCol()].deref_key(0);
        coor[1] = ite.getRow();
        row_nbmax_hcol(coor)++;

        coor[0] = _target.rowhid.data[ite.getRow()].deref_key(0);
        coor[1] = ite.getCol();
        col_nbmax_hrow(coor)++;

        coor[1] = _target.colhid.data[ite.getCol()].deref_key(0);
        maxgrid_nonzero_count(coor)++;

    }while (ite++);

    }

    uint32_t operator()(uint32_t thrID){
    ExOp::toZero(parascope[thrID].ll);
    parascope[thrID].dd_massive.toZero();
    parascope[thrID].dll_dMR2.toZero();
    parascope[thrID].nonzero_diff.toZero();
    parascope[thrID].nbmax_diff.toZero();
    parascope[thrID].indexes_changes.toMemfree();
    Tuple<uint32_t,2u> coor, acoor;

    uint32_t i,j;
    Tuple<Tuple<double,2u> > curproj[4]; // R H_j C_j
    double tmp;
    Tuple< Tuple<double, 2u>, 0u> future_statistics[3];
    Trianglix< Tuple<double, 3u>, 0u> future_statistic_trix[3];
    Tuple<double, 3u> trix_input;
    SparseTuple< Tuple<double, 2u> > future_factor;

    SparseTuple<double> candidates[2];
    SparseTuple<double> can_backup;
    Tuple<double, 3u> lls;
    Tuple<double, 3u> lld;
    Tuple<double, 2u> coeff;
    Tuple<double, 2u> coeff_test;
    double sum[4];

    FunctionScope<double&,double&,double> funabsoft(ExCo<double>::toAbsoftInverse,pow(0.5, 256.0f));

    Tuple<double, 2u> d_scale, o_scale, d_scale_tmp;
    Trianglix<double, 2u> dd_scale;

    Trianglix<double, 2u> dd_scale_tmp;

    Tuple<Tuple<double,2u>, 0u> rever;

    Tuple<double> prior_deriv;
    double prior_llterm[3];
    double Rval,Mval;

    double logexpmin;
    double expRterm;
    union LocalDef{
        double evilterms[6];
        KeyElem<double, Tuple<double, 5u> > evilterm_alt;
        LocalDef(){}
    };
    LocalDef ld;

    target.par_RM.getDims(coor);
    Tuple<Tuple<double,2u> > RM_deriv;

    Tuple<Tuple<Tuple<double,2u> > > rever_buf;
    Tuple<Tuple<double, 6u> > dd_buf;

    uint32_t dir;
    bool haschange;

    int insist=0;

    double prediction;

    double accuracy =0.0;
    parascope[thrID].nb_fail =0u;
    parascope[thrID].nb_insist =0u;
    uint32_t nb_funky=0u;

    Tuple<double, 4u> tosolve_tmp[3];

    Tuple<double, 3u> tosolve_num2;
    Tuple<double, 3u> tosolve_num2_tmp;
    Trianglix<double, 3u> tosolve_den2;
    Trianglix<double, 3u> tosolve_den2_tmp;
    Tuple<double> drop_buf;

    double expexp;
    double expexp2;
    bool checksafe;
    Vector<int> torem;
    uint32_t toremdirection;
    ProgressBarPrint pbar(20);
    if (parascope[thrID].use_stdout) pbar.start((istimefor_columns) ? "Compute columns": "Compute rows");

    // x_ij = NB ( hRz , hMz )  + NB(SoupFact_j * F_i, M_i)

    //

    if (istimefor_columns){

        RM_deriv.setSize(target.par_RM.sizes[1]);
        prior_deriv.setSize(target.par_RM.sizes[1]);

        rever_buf.setSize(target.rowhid.getNBcols());
        dd_buf.setSize(target.rowhid.getNBcols());


        for(i=0;i<3;i++) {
            future_statistics[i].setSize(target.par_RM.sizes[0]);
            future_statistic_trix[i].setSize(target.par_RM.sizes[0]);
        }
        ExOp::toZero(sum);
        for(coor[1]=colrange[thrID];coor[1] < colrange[thrID+1];coor[1]++){
            if (parascope[thrID].use_stdout) pbar.update(((double)coor[1] - colrange[thrID])/ (colrange[thrID+1]-colrange[thrID]));
            //printf("cur coor %i/%i\n", coor[1], target.colhid.getNBcols());fflush(stdout);
            const SparseTuple<double>& curcol = target.colhid.getColumn(coor[1]);
            if (curcol.getSize() == 0) {printf("Queried collumn %i is empty... impossible!\n",coor[1]); exit(1);}
            if (!ExOp::isValid(target.par_colscale[coor[1]][0])){
                newscale[coor[1]] = target.par_colscale[coor[1]];

                if (data.getColSize(coor[1]) == 0) {candidates[0] = target.colhid.getColumn(coor[1]); newhidden.memmoveColumn(coor[1], candidates[0]); continue;} // nothing to see...
                // column only has 1s, cannot be used to update R or M, drop out only, pure state guarrantied

                // put in stone... for now
                candidates[0] = target.colhid.getColumn(coor[1]); newhidden.memmoveColumn(coor[1], candidates[0]); continue;
                // put in stone... for now

                RM_deriv.toZero();

                for(coor[0]=0;coor[0]<col_nbmax_hrow.dims[0];coor[0]++){
                    acoor[0] = coor[0];
                    for(acoor[1] =0;acoor[1]< target.par_RM.sizes[1]; acoor[1]++){
                        RM_deriv[acoor[1]][0] += col_nbmax_hrow(coor) * log(target.par_D(acoor));
                        RM_deriv[acoor[1]][0] += (nb_maxrowh[coor[0]] - col_nbmax_hrow(coor)) * log(1.0 - target.par_D(acoor));
                    }
                }
                acoor[1] =0;
                for(acoor[0]=1 ;acoor[0]< target.par_RM.sizes[1]; acoor[0]++){
                    if (RM_deriv[acoor[0]][0] > RM_deriv[acoor[1]][0]) acoor[1] = acoor[0];
                }
                candidates[0].toZero();
                candidates[0][acoor[1]] = 1.0;
                if (candidates[0].deref_key(0) != curcol.deref_key(0)){
                    parascope[thrID].nbmax_diff[candidates[0].deref_key(0)]++;
                    parascope[thrID].nbmax_diff[curcol.deref_key(0)]--;
                    for(i = 0; i < data.data[coor[1]].getSize();i++){
                        acoor[0] = target.rowhid.data[ data.data[coor[1]].deref_key(i)].deref_key(0);
                        acoor[1] = curcol.deref_key(0);
                        parascope[thrID].nonzero_diff(acoor)--;
                        acoor[1] = candidates[0].deref_key(0);
                        parascope[thrID].nonzero_diff(acoor)++;
                    }
                }
                newhidden.memmoveColumn(coor[1], candidates[0]);
            }else{
                checksafe = false;

                prior_deriv = target.colhid_prior * curcol;

                if (coor[1] >= target.colhid.data.getSize()) exit(1);
                // Step 1 cache line specific data
                curproj[0] = target.par_RM * curcol;

                //dll_dH.toZero();
                //dll_dC.toZero();
                // Step 2 find direction, find 2 vector for R and M
                RM_deriv.toZero();



                prior_llterm[0] = target.colhid_prior.mkInnerMult(curcol)+ target.colhid_prior_constant;
                lls[2] = log(prior_llterm[0]) * ((tmptmp) ? data.data[coor[1]].getSize() : 1);
                //lls[2] = 0.0;
                future_statistics[2].toZero(); future_statistic_trix[2].toZero();
                d_scale.toZero(); dd_scale.toZero();
                for(i=0;i<data.data[coor[1]].getSize();i++){
                    coor[0]= data.data[coor[1]].deref_key(i);
                    if (!ExOp::isValid(target.par_rowscale[coor[0]])) continue;
                    coeff = curproj[0].mkInnerProd(target.rowhid.getColumn(coor[0]));
                    coeff += target.par_rowscale[coor[0]];
                    coeff += target.par_colscale[coor[1]];
                    deflatedNBterm_wr_f_d2(ld.evilterms, data.data[coor[1]].deref(i)  , coeff[0], coeff[1]);
                    dd_scale.data[0] += ld.evilterms[3]; dd_scale.data[1] += ld.evilterms[5]; dd_scale.data[2] += ld.evilterms[4];
                    d_scale[0] += ld.evilterms[1]; d_scale[1] += ld.evilterms[2];

                    lls[2] += ld.evilterms[0];
                    //printf("(%e,%e), %e,%e,%e,%e,%e\n", target.par_rowscale[coor[0]][0], target.par_rowscale[coor[0]][1], ld.evilterms[0],ld.evilterms[1],ld.evilterms[2],ld.evilterms[3],ld.evilterms[4] );

                    //d_scale[0] += e += ld.evilterm_alt;
                    rever = target.par_RM.mkBackMult(target.rowhid.getColumn(coor[0]));

                    if (auto ite = rever.getIterator()) do{
                        RM_deriv[ite()][0] += (*ite)[0] * ld.evilterms[1];
                        RM_deriv[ite()][1] += (*ite)[1] * ld.evilterms[2];
                    }while(ite++);

                    if (auto ite = target.rowhid.getColumn(coor[0]).getIterator()) do{

                        future_statistics[2][ite()][0] += ld.evilterm_alt.d[0] * (*ite);
                        future_statistics[2][ite()][1] += ld.evilterm_alt.d[1] * (*ite);
                        trix_input[0] = ld.evilterms[3] * (*ite);
                        trix_input[1] = ld.evilterms[4] * (*ite);
                        trix_input[2] = ld.evilterms[5] * (*ite) * 0.5;

                        if (auto ite2 = target.rowhid.getColumn(coor[0]).getIterator()) do{
                            if (ite2() > ite()) continue;
                            future_statistic_trix[2].data[ite2() + ((ite() * (ite() + 1)) >> 1)] += trix_input * (*ite2);
                        }while(ite2++);
                    }while(ite++);

                    rever_buf[i].toMemmove(rever);
                    dd_buf[i][0] = ld.evilterms[0];
                    dd_buf[i][1] = ld.evilterms[1];
                    dd_buf[i][2] = ld.evilterms[2];
                    dd_buf[i][3] = ld.evilterms[3];
                    dd_buf[i][4] = ld.evilterms[4];
                    dd_buf[i][5] = ld.evilterms[5];

                } // dG/dz = c1f1 VM  dH/dz = c2 f2VR
                //printf("cur coor %i/%i\n", coor[1], target.colhid.getNBcols());fflush(stdout);

                parascope[thrID].ll[0] += parascope[thrID].lls[2];
                prediction = 0.0f;
                // Step 3 (line search) define set of candidates

                sum[0] = 0.0;
                prior_llterm[2] = 0.0;
                candidates[0].toMemfree();
                for(i=0;i< RM_deriv.getSize();i++){
                    prior_llterm[1] = RM_deriv[i][0] + RM_deriv[i][1] + prior_deriv[i] * (2.0 * ((tmptmp) ? data.data[coor[1]].getSize() : 1) / prior_llterm[0]);
                    if ((curcol.find(i) != 0xFFFFFFFF)||(prior_llterm[1] > 0.0f)) {
                        candidates[0][i] = prior_llterm[1];
                        sum[0] +=  prior_llterm[1];
                    }
                }
                if (candidates[0].getSize() < 2u) { // stuck in a corner,... that's good but...
                    candidates[0] = curcol;
                    dd_scale_tmp = dd_scale; d_scale_tmp = d_scale;
                    d_scale = dd_scale.toEigenTransform(funabsoft) * d_scale;
                    d_scale *= hiddenalpha;
                }else{
                    // project onto plane
                    coor[0]=0;
                    do{
                        torem.toMemfree();
                        sum[0] /= candidates[0].getSize();
                        if (auto ite = candidates[0].getIterator()) do{
                            if (((*ite) > sum[0])||(curcol.find(ite()) != 0xFFFFFFFF)) (*ite) -= sum[0];
                            else torem.push_back(ite());
                        }while(ite++);
                        if (torem.getSize() == 0u) {
                            sum[0] = 0.0;
                            if (auto ite = candidates[0].getIterator()) do{sum[0] += (*ite);}while(ite++);
                            if (fabs(sum[0]) < 0.0001) break;
                            if ((coor[0]++) == 5) {
                                coor[0] = 0;
                                if (auto ite = candidates[0].getIterator()) do{if ((*ite)>= 0) coor[0]++;}while(ite++);
                                if (auto ite = candidates[0].getIterator()) do{(*ite) = ((*ite)>= 0) ?  candidates[0].getSize() - coor[0]: coor[0];}while(ite++);
                                if (auto ite = candidates[0].getIterator()) do{(*ite) /=  candidates[0].getSize();}while(ite++);
                                //printf("warning: got stuck... catastrophy cancellations\n");
                                break;
                            }
                        }else{
                            if (candidates[0].getSize() - torem.getSize() <= 1){ // would leave nothing... try the highest second?
                                if (candidates[0].getSize() == torem.getSize()) {printf("Impossible!\n"); exit(1);}
                                acoor[0] = 0;
                                for(acoor[1]=1;acoor[1] < torem.getSize();acoor[1]++) if (candidates[0][torem[acoor[0]]] < candidates[0][torem[acoor[1]]]) acoor[0] = acoor[1];
                                for(i=0;i<torem.getSize();i++) {
                                    if (i == acoor[0]) continue;
                                    j = candidates[0].find(torem[i]);
                                    if (j == 0xFFFFFFFF) {printf("erasing unexisting... !\n"); exit(1);}
                                    candidates[0].erase_from_iterator(j);
                                }
                                 if (candidates[0].getSize() != 2) {printf("2 should remain. got %i\n", candidates[0].getSize()); exit(1);}
                                if (candidates[0].deref_key(0) == torem[acoor[0]]){
                                    candidates[0].deref(0) = 1.0;
                                    candidates[0].deref(1) = -1.0;
                                }else{
                                    candidates[0].deref(0) = -1.0;
                                    candidates[0].deref(1) = 1.0;
                                }
                                break;
                            }
                            for(i=0;i<torem.getSize();i++) {
                                j = candidates[0].find(torem[i]);
                                if (j == 0xFFFFFFFF) {printf("trying to remove un-existing!\n"); exit(1);}
                                candidates[0].erase_from_iterator(j);
                            }
                            sum[0] =0.0;
                            if (auto ite = candidates[0].getIterator()) do{sum[0] += *ite;}while(ite++);
                        }
                    }while(true);

                    if (auto ite = candidates[0].getIterator()){
                        while((*ite) > 0.0) ite++;
                        sum[0] = curcol[dir = ite()] / (*ite);

                        while(ite++){
                            if (*ite > 0.0) continue;
                            sum[1] = curcol[ite()] / (*ite);
                            if (sum[1] > sum[0]) {sum[0] = sum[1]; dir = ite();}
                        }
                        sum[0] = -sum[0]; // w has to be in 0-sum[0] range

                        tosolve_num2[0] = (target.colhid_prior.mkInnerMult(candidates[0],curcol)+ target.colhid_prior_constant) * ( 2.0 * ((tmptmp) ? data.data[coor[1]].getSize() : 1)/ prior_llterm[0]);
                        tosolve_num2[1] = d_scale[0];
                        tosolve_num2[2] = d_scale[1];
                        tosolve_den2.data[1] = 0.0; tosolve_den2.data[3] = 0.0;
                        tosolve_den2.data[2] = dd_scale.data[0];
                        tosolve_den2.data[5] = dd_scale.data[2];
                        tosolve_den2.data[4] = dd_scale.data[1];
                        tosolve_den2.data[0] = (target.colhid_prior.mkInnerMult(candidates[0])+ target.colhid_prior_constant) * ( 2.0 * ((tmptmp) ? data.data[coor[1]].getSize() : 1) / prior_llterm[0]);
                        for(i=0;i<data.data[coor[1]].getSize();i++){
                            coor[0]= data.data[coor[1]].deref_key(i);
                            if (!ExOp::isValid(target.par_rowscale[coor[0]])) continue;
                            coeff = rever_buf[i].mkInnerProd(candidates[0]);
                            tosolve_num2[0] += coeff[0] * dd_buf[i][1]+ coeff[1] * dd_buf[i][2];
                            tosolve_den2.data[0] += dd_buf[i][3] * coeff[0] * coeff[0] + dd_buf[i][4] * coeff[1] * coeff[1] + dd_buf[i][5] * coeff[0] * coeff[1];
                            tosolve_den2.data[1] += dd_buf[i][3] * coeff[0] + dd_buf[i][5] * coeff[1] * 0.5;
                            tosolve_den2.data[3] += dd_buf[i][4] * coeff[1] + dd_buf[i][5] * coeff[0] * 0.5;
                        }

                        tosolve_den2_tmp = tosolve_den2;
                        tosolve_num2_tmp = tosolve_num2;
                        tosolve_den2.toEigenTransform(funabsoft);
                        tosolve_num2 = tosolve_den2 * tosolve_num2;

                        tosolve_num2 *= hiddenalpha;

                        expexp = tosolve_num2[0] * tosolve_num2_tmp[0] + tosolve_num2[1] * tosolve_num2_tmp[1] + tosolve_num2[2] * tosolve_num2_tmp[2];
                        expexp2 = tosolve_den2_tmp.data[0] * tosolve_num2[0] * tosolve_num2[0] * 0.5;
                        expexp2 += tosolve_den2_tmp.data[2] * tosolve_num2[1] * tosolve_num2[1] * 0.5;
                        expexp2 += tosolve_den2_tmp.data[5] * tosolve_num2[2] * tosolve_num2[2] * 0.5;
                        expexp2 += tosolve_den2_tmp.data[1] * tosolve_num2[0] * tosolve_num2[1];
                        expexp2 += tosolve_den2_tmp.data[3] * tosolve_num2[0] * tosolve_num2[2];
                        expexp2 += tosolve_den2_tmp.data[4] * tosolve_num2[1] * tosolve_num2[2];

                        if (tosolve_num2[0] <= 0.0){ // went the other way... just update scale
                            if (tosolve_num2[0] == 0.0){
                                candidates[0] = curcol;
                            }else{
                                candidates[0] = curcol;
                                dd_scale_tmp = dd_scale; d_scale_tmp = d_scale;
                                d_scale = dd_scale.toEigenTransform(funabsoft) * d_scale;
                                d_scale *= hiddenalpha;

                            }
                       }else if (tosolve_num2[0] >= sum[0]){ // too far!, find solution on boundary...
                            d_scale[0] = tosolve_num2_tmp[1] + sum[0] * tosolve_den2.data[1];
                            d_scale[1] = tosolve_num2_tmp[2] + sum[0] * tosolve_den2.data[3];
                            d_scale = dd_scale.toEigenTransform(funabsoft) * d_scale;
                            d_scale *= hiddenalpha;

                            candidates[1] = candidates[0];
                            if (auto ite = candidates[0].getIterator()) do{
                                (*ite) = (*ite) * sum[0]  + (((j = curcol.find(ite())) == 0xFFFFFFFF) ? 0.0: curcol.deref(j));
                            }while(ite++);

                            if (auto ite = curcol.getIterator()) do{
                                if (candidates[0].find(ite())== 0xFFFFFFFF) candidates[0][ite()] = *ite;
                            }while(ite++);
                            j = candidates[0].find(dir);
                            if (j == 0xFFFFFFFF) {printf("trying to remove un-existing!\n"); exit(1);}
                            candidates[0].erase_from_iterator(j);
                            sum[1] =0.0;
                            if (auto ite = candidates[0].getIterator()) {
                                do{
                                    sum[1] += *ite;
                                }while(ite++);
                                if (fabs(sum[1] - 1.0) > 0.1) {printf("bound proj dir %i... hidden state norm is %e! %e is fact\n", dir, sum[1], sum[0]);
                                candidates[1].show();
                                candidates[0].show();
                                curcol.show();
                                LFH_ALIVE;exit(1);}
                            }else {printf("got empty candidate!\n"); exit(1);}
                            candidates[1].toMemfree();
                        }else{ // all set!
                            //checksafe = true;

                            d_scale[0] = tosolve_num2[1];
                            d_scale[1] = tosolve_num2[2];
                            candidates[1] = candidates[0];
                            if (auto ite = candidates[0].getIterator()) do{
                                (*ite) = (*ite) * tosolve_num2[0] + (((j = curcol.find(ite())) == 0xFFFFFFFF) ? 0.0: curcol.deref(j));
                            }while(ite++);
                            if (auto ite = curcol.getIterator()) do{
                                if (candidates[0].find(ite())== 0xFFFFFFFF) candidates[0][ite()] = *ite;
                            }while(ite++);
                            sum[1] = 0.0;
                            if (auto ite = candidates[0].getIterator()) {
                                do{
                                    sum[1] += *ite;
                                }while(ite++);
                                if (fabs(sum[1] - 1.0) > 0.1) {printf("normal proj... hidden state norm is %e! %e is fact\n", sum[1], tosolve_num2[0]);
                                candidates[0].show();
                                candidates[1].show();
                                curcol.show();
                                LFH_ALIVE; exit(1);}
                            }else {printf("got empty candidate!\n"); exit(1);}
                            candidates[1].toMemfree();
                        }
                        candidates[0].makeMaximumAsFirst();
                        if (!ignore_drop){
                            if (candidates[0].deref_key(0) != curcol.deref_key(0)) {
                                // major change, add *last* boundary on the search space!
                                tosolve_num2[2] = ((i = curcol.find(candidates[0].deref_key(0))) == 0xFFFFFFFF) ? 0.0 : curcol.deref(i);
                                if (auto ite = curcol.getIterator()) {
                                    tosolve_num2[0] = (*ite) - tosolve_num2[2]; // guarranied to be positive
                                    tosolve_num2[1] = candidates[0].deref(0) -  (((i = candidates[0].find(ite())) == 0xFFFFFFFF) ? 0.0 : candidates[0].deref(i));
                                    sum[1] = tosolve_num2[1] / (tosolve_num2[0] + tosolve_num2[1]);
                                    j = 0;
                                    while(ite++){
                                        if ((tosolve_num2[0] =  (*ite) - tosolve_num2[2]) <= 0) continue;
                                        tosolve_num2[1] = candidates[0].deref(0) -  (((i = candidates[0].find(ite())) == 0xFFFFFFFF) ? 0.0 : candidates[0].deref(i));
                                        sum[0] = tosolve_num2[1] / (tosolve_num2[0] + tosolve_num2[1]);
                                        if (sum[0] < sum[1]) {j = ite.i; sum[1] = sum[0];}
                                    }
                                }
                                sum[0] = candidates[0].deref(0) * (1.0 - sum[1]) + sum[1] * (((i = curcol.find(candidates[0].deref_key(0))) == 0xFFFFFFFF) ? 0.0 : curcol.deref(i));
                                candidates[1][curcol.deref_key(j)] = sum[0];
                                candidates[1][candidates[0].deref_key(0)] = sum[0];
                                if (auto ite = curcol.getIterator()) do {
                                    if (candidates[1].find(ite()) == 0xFFFFFFFF) candidates[1][ite()] = (*ite) * sum[1] + (1.0 - sum[1]) * (((i = candidates[0].find(ite())) == 0xFFFFFFFF) ? 0.0 : candidates[0].deref(i));
                                }while(ite++);
                                if (auto ite = candidates[0].getIterator()) do {
                                    if (candidates[1].find(ite()) == 0xFFFFFFFF) candidates[1][ite()] = (*ite) * (1.0 - sum[1]);
                                }while(ite++);

                                if (auto ite = candidates[1].getIterator()) {
                                    sum[0] = 0.0;
                                    do{
                                        sum[0] += *ite;
                                    }while(ite++);
                                    if (fabs(sum[0] - 1.0) > 0.1) {printf("hidded state norm cand1 is %e!\n", sum[0]);
                                        candidates[1].show();
                                        printf("frac %e\n", sum[1]);
                                        candidates[0].show();
                                        curcol.show();
                                    LFH_ALIVE;exit(1);}
                                }
                            }/*else{ // forcefully try closest maybe?

                            }*/
                        }
                    }else{
                        // stuck in corner again...
                        candidates[0] = curcol;
                        dd_scale_tmp = dd_scale; d_scale_tmp = d_scale;
                        d_scale = dd_scale.toEigenTransform(funabsoft) * d_scale;
                        d_scale *= hiddenalpha;
                    }
                }
                sum[0] = 0.0f;
                if (auto ite = candidates[0].getIterator()) {
                    do{
                        sum[0] += *ite;
                    }while(ite++);
                    if (fabs(sum[0] - 1.0) > 0.1) {printf("hidded state norm is hh %e!\n", sum[0]);
                        candidates[0].show();
                        candidates[1].show();
                    exit(1);}
                }else {printf("got empty candidate!\n"); exit(1);}

                // find intersection on line
                // new scale!
                //printf("cur coor %i/%i\n", coor[1], target.colhid.getNBcols());fflush(stdout);

                coeff = target.par_colscale[coor[1]];
                coeff += d_scale;
                if ((tmp = (coeff[0]*coeff[0] + coeff[1]*coeff[1])) > 100.0){ // too far, send to boundary
                    tmp = pow(tmp, -0.5);
                    d_scale[0] = coeff[0] * tmp;
                    d_scale[1] = coeff[1] * tmp;
                    d_scale -= target.par_colscale[coor[1]];
                }

                // Step 4 Evaluate Quantidates

                insist =0;
                while(true){

                    for(j=0;j<2;j++){
                        future_statistics[j].toZero(); future_statistic_trix[j].toZero();
                        curproj[j] = target.par_RM * candidates[j];
                        lls[j] = log(target.colhid_prior.mkInnerMult(candidates[j])+ target.colhid_prior_constant) * ((tmptmp) ? data.data[coor[1]].getSize() : 1);
                       // lls[j] = 0.0;
                        if ((candidates[1].getSize() == 0)||(insist > 0)) break;
                    }
                    for(i=0;i<data.data[coor[1]].getSize();i++){
                        coor[0]= data.data[coor[1]].deref_key(i);
                        if (!ExOp::isValid(target.par_rowscale[coor[0]])) continue;
                        for(j=0;j<2;j++){
                            coeff = curproj[j].mkInnerProd(target.rowhid.getColumn(coor[0]));
                            coeff += target.par_rowscale[coor[0]];
                            coeff += target.par_colscale[coor[1]];
                            coeff += d_scale;


                            deflatedNBterm_wr_f_d2(ld.evilterms, data.data[coor[1]].deref(i), coeff[0], coeff[1]);
                            if (auto ite = target.rowhid.getColumn(coor[0]).getIterator()) do{
                                future_statistics[j][ite()][0] += ld.evilterm_alt.d[0] * (*ite);
                                future_statistics[j][ite()][1] += ld.evilterm_alt.d[1] * (*ite);
                                trix_input[0] = ld.evilterms[3] * (*ite);
                                trix_input[1] = ld.evilterms[4] * (*ite);
                                trix_input[2] = ld.evilterms[5] * (*ite) * 0.5;
                                if (auto ite2 = target.rowhid.getColumn(coor[0]).getIterator()) do{
                                    if (ite2() > ite()) continue;
                                    future_statistic_trix[j].data[ite2() + ((ite() * (ite() + 1)) >> 1)] += trix_input * (*ite2);
                                }while(ite2++);
                            }while(ite++);

                            lls[j] += ld.evilterms[0];
                            if ((candidates[1].getSize() == 0)||(insist > 0)) break;
                        }
                    } // H(z) = fcfc vMMv

                    if (checksafe){
                        double diffdiff = (log(target.colhid_prior.mkInnerMult(candidates[0])+ target.colhid_prior_constant) - log(target.colhid_prior.mkInnerMult(curcol)+ target.colhid_prior_constant)) * ((tmptmp) ? data.data[coor[1]].getSize() : 1);
                        printf("%e vs %e predicted\n", lls[0] - lls[2], expexp + expexp2);
                        printf("quad term: %e vs %e\n", lls[0] - lls[2] - expexp, expexp2);

                    }

                    // Step 5 select best or do not update
                    //tmp = d_scale_tmp[0] * d_scale[0] + d_scale_tmp[1] * d_scale[1];
                    //tmp += (dd_scale_tmp.data[0] * d_scale[0] * d_scale[0] + dd_scale_tmp.data[2] * d_scale[1] * d_scale[1])* 0.5 + dd_scale_tmp.data[1] * d_scale[0] * d_scale[1];
                    //printf("%e incr (vs%e)\n", lls[0] - lls[3], tmp);

                //    for(coor[0]=0;coor[0] < new_nbmax_h.getSize();coor[0]++) new_nbmax_h(coor) = row_nbmax_hcol(coor);
                    if (candidates[1].getSize() == 0){
                         j = (lls[0] > lls[2]) ? 0 : 2;
                    }else{ // lls[0], lls[1] lls[3]
                        // worrying about dropouts!
                        lld.toZero();
                        acoor[1] = curcol.deref_key(0);
                        for(acoor[0] =0;acoor[0]< target.par_RM.sizes[0]; acoor[0]++){coor[0] = acoor[0];
                            lld[2] += col_nbmax_hrow(coor) * log(target.par_D(acoor));
                            lld[2] += (nb_maxrowh[coor[0]] - col_nbmax_hrow(coor)) * log(1.0 - target.par_D(acoor));
                        }
                        acoor[0] = candidates[0].deref_key(0);
                        for(acoor[0] =0;acoor[0]< target.par_RM.sizes[0]; acoor[0]++){coor[0] = acoor[0];
                            lld[0] += col_nbmax_hrow(coor) * log(target.par_D(acoor));
                            lld[0] += (nb_maxrowh[coor[0]] - col_nbmax_hrow(coor)) * log(1.0 - target.par_D(acoor));
                        }
                        if (candidates[1].deref_key(0) == curcol.deref_key(0)) lld[1] = lld[2];
                        else{
                            acoor[0] = candidates[1].deref_key(0);
                            for(acoor[0] =0;acoor[0]< target.par_RM.sizes[0]; acoor[0]++){coor[0] = acoor[0];
                                lld[1] += col_nbmax_hrow(coor) * log(target.par_D(acoor));
                                lld[1] += (nb_maxrowh[coor[0]] - col_nbmax_hrow(coor)) * log(1.0 - target.par_D(acoor));
                            }
                        }

                        j = (lls[0] + lld[0] > lls[2] + lld[2]) ? 0 : 2;
                        if ((lls[1] + ((lld[0] > lld[1]) ? lld[0]  : lld[1])) > lls[j] + lld[j]) {
                            if (lld[0] > lld[1]) candidates[1].swapEntries(0,1);
                            j=1;
                        }
                        if ((j != 2)&&(candidates[j].deref_key(0) != curcol.deref_key(0))){
                            // registering change...
                            parascope[thrID].nbmax_diff[candidates[j].deref_key(0)]++;
                            parascope[thrID].nbmax_diff[curcol.deref_key(0)]--;
                            for(i = 0; i < data.data[coor[1]].getSize();i++){
                                acoor[0] = target.rowhid.data[ data.data[coor[1]].deref_key(i)].deref_key(0);
                                acoor[1] = curcol.deref_key(0);
                                parascope[thrID].nonzero_diff(acoor)--;
                                acoor[1] = candidates[j].deref_key(0);
                                parascope[thrID].nonzero_diff(acoor)++;
                            }
                        }
                    }

                    if (j != 2) break;
                    if (insist++ > 10) break;
                    if (auto ite = curcol.getIterator()) do{ candidates[0][ite()] += *ite; }while(ite++);
                    if (auto ite = candidates[0].getIterator()) do{ *ite *= 0.5; }while(ite++);
                    d_scale *= 0.5;

                    candidates[0].makeMaximumAsFirst();
                    if (candidates[0].deref_key(0) == curcol.deref_key(0)) candidates[1].toMemfree();

                } // insist loop end
                parascope[thrID].nb_insist += parascope[thrID].insist;

                // Step 6 save scope for parameter update
                if (j == 2) { // failed to improve...
                    //printf("col %i failed  %e %e\n", coor[1], target.par_colscale[coor[1]][0], target.par_colscale[coor[1]][1] );
                    candidates[0] = curcol;
                    newhidden.memmoveColumn(coor[1], candidates[0]);
                    newscale[coor[1]] = target.par_colscale[coor[1]];
                    parascope[thrID].nb_fail++;
                }else{
                    newhidden.memmoveColumn(coor[1], candidates[j]);
                    newscale[coor[1]] = target.par_colscale[coor[1]] + d_scale; // to change?
                }
                parascope[thrID].ll[1] += parascope[thrID].lls[j];
                candidates[1].toMemfree();

                future_factor.toMemfree();
                if (auto ite = newhidden.data[coor[1]].getIterator()) do{
                    future_factor[ite()][0] = *ite; future_factor[ite()][1] = *ite;
                }while(ite++);

                parascope[thrID].dll_dMR2.toAddOuterProd(future_statistics[j], future_factor);

                uint32_t dashift = target.par_RM.sizes[0] * target.par_RM.sizes[1];
                if (auto ite = newhidden.data[coor[1]].getIterator()) do{
                    if (auto ite2 = newhidden.data[coor[1]].getIterator()) do{
                        if (ite2() > ite()) continue;
                        tmp = (*ite2) * (*ite);
                        uint32_t k;
                        for(i = 0,k=0; i< target.par_RM.sizes[0];i++){
                            uint32_t l = i + ite() * target.par_RM.sizes[0];

                            double* currow_A = parascope[thrID].dd_massive.data + (((l * (l +1)) >> 1) + ite2() * target.par_RM.sizes[0]);
                            double* currow_B = parascope[thrID].dd_massive.data + ((((l + dashift) * (l + dashift +1)) >> 1) + ite2() * target.par_RM.sizes[0]);
                            for(l=0;l<=i;l++,k++){
                                currow_A[l] += future_statistic_trix[j].data[k][0] * tmp;
                                currow_B[l] += future_statistic_trix[j].data[k][2] * tmp;
                                currow_B[l + dashift] += future_statistic_trix[j].data[k][1] * tmp;
                            }
                        }
                        //printf("%i is k  %i\n", k,  future_statistic_trix[j].getSize());
                    }while(ite2++);
                }while(ite++);
            }
        }
        lls[0] =0;
        for(coor[1]=rowrange[thrID];coor[1] < rowrange[thrID+1];coor[1]++){
            lls[0] += log(target.rowhid_prior.mkInnerMult(target.rowhid.getColumn(coor[1]))+ target.rowhid_prior_constant) * ((tmptmp) ? trdata.data[coor[1]].getSize() : 1);
        }
        parascope[thrID].ll[0] += lls[0];
        parascope[thrID].ll[1] += lls[0];

    }else{ // time for rows!
        RM_deriv.setSize(target.par_RM.sizes[0]);
        prior_deriv.setSize(target.par_RM.sizes[0]);


        rever_buf.setSize(target.colhid.getNBcols());
        dd_buf.setSize(target.colhid.getNBcols());


        for(i=0;i<3;i++) {
            future_statistics[i].setSize(target.par_RM.sizes[1]);
            future_statistic_trix[i].setSize(target.par_RM.sizes[1]);
        }
        ExOp::toZero(sum);

        for(coor[1]=rowrange[thrID];coor[1] < rowrange[thrID+1];coor[1]++){
            if (parascope[thrID].use_stdout) pbar.update(((double)coor[1]-rowrange[thrID])/ (rowrange[thrID+1]-rowrange[thrID]));
            if (trdata.data[coor[1]].getSize() > target.colhid.getNBcols()) {printf("Holly Molly! %i and %i",trdata.data[coor[1]].getSize(), target.colhid.getNBcols());exit(1);}
            //printf("cur coor %i/%i\n", coor[1], target.colhid.getNBcols());fflush(stdout);
            const SparseTuple<double>& curcol = target.rowhid.getColumn(coor[1]);
            if (curcol.getSize() == 0) {printf("Queried collumn %i is empty... impossible!\n",coor[1]); exit(1);}
            if (!ExOp::isValid(target.par_rowscale[coor[1]][0])){
                newscale[coor[1]] = target.par_rowscale[coor[1]];

                if (trdata.getColSize(coor[1]) == 0) {candidates[0] = target.rowhid.getColumn(coor[1]);newhidden.memmoveColumn(coor[1], candidates[0]); continue;} // nothing to see...
                // column only has 1s, cannot be used to update R or M, drop out only, pure state guarrantied

                // put in stone... for now
                candidates[0] = target.rowhid.getColumn(coor[1]);newhidden.memmoveColumn(coor[1], candidates[0]);
                // put in stone... for now

                RM_deriv.toZero();

                for(coor[0]=0;coor[0]<row_nbmax_hcol.dims[0];coor[0]++){acoor[1] = coor[0];
                    for(acoor[0] =0;acoor[0]< target.par_RM.sizes[0]; acoor[0]++){
                        RM_deriv[acoor[0]][0] += row_nbmax_hcol(coor) * log(target.par_D(acoor));
                        RM_deriv[acoor[0]][0] += (nb_maxcolh[coor[0]] - row_nbmax_hcol(coor)) * log(1.0 - target.par_D(acoor));
                    }
                }
                acoor[1] =0;
                for(acoor[0]=0 ;acoor[0]< target.par_RM.sizes[0]; acoor[0]++){
                    if (RM_deriv[acoor[0]][0] > RM_deriv[acoor[1]][0]) acoor[1] = acoor[0];
                }
                candidates[0].toZero();
                candidates[0][acoor[1]] = 1.0;
                if (candidates[0].deref_key(0) != curcol.deref_key(0)){
                    parascope[thrID].nbmax_diff[candidates[0].deref_key(0)]++;
                    parascope[thrID].nbmax_diff[curcol.deref_key(0)]--;
                    for(i = 0; i < trdata.data[coor[1]].getSize();i++){
                        acoor[1] = target.colhid.data[trdata.data[coor[1]].deref_key(i)].deref_key(0);
                        acoor[0] = curcol.deref_key(0);
                        parascope[thrID].nonzero_diff(acoor)--;
                        acoor[0] = candidates[0].deref_key(0);
                        parascope[thrID].nonzero_diff(acoor)++;
                    }
                }
                newhidden.memmoveColumn(coor[1], candidates[0]);
            }else{

                prior_deriv = target.rowhid_prior * curcol;
                // Step 1 cache line specific data
                curproj[0] = target.par_RM.mkBackMult(curcol);

                //dll_dH.toZero();
                //dll_dC.toZero();

                // Step 2 find direction, find 2 vector for R and M
                RM_deriv.toZero();
                prior_llterm[0] = target.rowhid_prior.mkInnerMult(curcol) + target.rowhid_prior_constant;
                lls[2] = log(prior_llterm[0]) * ((tmptmp) ? trdata.data[coor[1]].getSize() : 1);

                future_statistics[2].toZero();future_statistic_trix[2].toZero();
                d_scale.toZero();dd_scale.toZero();
               // printf("loopya\n"); fflush(stdout);
                for(i=0;i<trdata.data[coor[1]].getSize();i++){
                    coor[0]= trdata.data[coor[1]].deref_key(i);
                    if (!ExOp::isValid(target.par_colscale[coor[0]][0])) continue; // can do better...
                    coeff = curproj[0].mkInnerProd(target.colhid.getColumn(coor[0]) );
                    coeff += target.par_colscale[coor[0]];
                    coeff += target.par_rowscale[coor[1]];
                    deflatedNBterm_wr_f_d2(ld.evilterms, trdata.data[coor[1]].deref(i), coeff[0], coeff[1]);
                    dd_scale.data[0] += ld.evilterms[3]; dd_scale.data[1] += ld.evilterms[5]; dd_scale.data[2] += ld.evilterms[4];
                    d_scale[0] += ld.evilterms[1]; d_scale[1] += ld.evilterms[2];
                    lls[2] += ld.evilterms[0];

                    rever = target.par_RM * target.colhid.getColumn(coor[0]);
                    if (auto ite = rever.getIterator()) do{
                        RM_deriv[ite()][0] += (*ite)[0] * ld.evilterms[1];
                        RM_deriv[ite()][1] += (*ite)[1] * ld.evilterms[2];
                    }while(ite++);


                    if (auto ite = target.colhid.getColumn(coor[0]).getIterator()) do{
                        future_statistics[2][ite()][0] += ld.evilterm_alt.d[0] * (*ite);
                        future_statistics[2][ite()][1] += ld.evilterm_alt.d[1] * (*ite);
                        trix_input[0] = ld.evilterms[3] * (*ite);
                        trix_input[1] = ld.evilterms[4] * (*ite);
                        trix_input[2] = ld.evilterms[5] * (*ite) * 0.5;
                        if (auto ite2 = target.colhid.getColumn(coor[0]).getIterator()) do{
                            if (ite2() > ite()) continue;
                            future_statistic_trix[2].data[ite2() + ((ite() * (ite() + 1)) >> 1)] += trix_input * (*ite2);
                        }while(ite2++);

                    }while(ite++);

                    rever_buf[i].toMemmove(rever);
                    dd_buf[i][0] = ld.evilterms[0];
                    dd_buf[i][1] = ld.evilterms[1];
                    dd_buf[i][2] = ld.evilterms[2];
                    dd_buf[i][3] = ld.evilterms[3];
                    dd_buf[i][4] = ld.evilterms[4];
                    dd_buf[i][5] = ld.evilterms[5];
                }
                //dd_scale.data[1] *= 0.5;

                parascope[thrID].ll[0] += lls[2];
                // Step 3 (line search) define set of candidates


                sum[0] = 0.0;
                prior_llterm[2] = 0.0;
                candidates[0].toZero();
                for(i=0;i< RM_deriv.getSize();i++){
                    prior_llterm[1] = RM_deriv[i][0] + RM_deriv[i][1] + prior_deriv[i] * (2.0 * ((tmptmp) ? trdata.data[coor[1]].getSize() : 1) / prior_llterm[0]);
                    if ((curcol.find(i) != 0xFFFFFFFF)||(prior_llterm[1] > 0.0f)) {
                        candidates[0][i] = prior_llterm[1];
                        sum[0] +=  prior_llterm[1];
                    }
                }
                if (candidates[0].getSize() < 2u) { // stuck in a corner,... that's good but...
                    candidates[0] = curcol;
                    dd_scale_tmp = dd_scale; d_scale_tmp = d_scale;
                    d_scale = dd_scale.toEigenTransform(funabsoft) * d_scale;
                    d_scale *= hiddenalpha;
                }else{
                    // project onto plane
                    coor[0] =0;

                    do{
                        torem.toMemfree();
                        sum[0] /= candidates[0].getSize();
                        if (auto ite = candidates[0].getIterator()) do{
                            if (((*ite) > sum[0])||(curcol.find(ite()) != 0xFFFFFFFF)) (*ite) -= sum[0];
                            else torem.push_back(ite());
                        }while(ite++);
                        if (torem.getSize() == 0u) {
                            sum[0] = 0.0;
                            if (auto ite = candidates[0].getIterator()) do{sum[0] += (*ite);}while(ite++);
                            if (fabs(sum[0]) < 0.0001) break;
                            if ((coor[0]++) == 5) {
                                printf("warning: got stuck... catastrophy cancellations\n");
                                coor[0] = 0;
                                if (auto ite = candidates[0].getIterator()) do{if ((*ite)>= 0) coor[0]++;}while(ite++);
                                if (auto ite = candidates[0].getIterator()) do{(*ite) = ((*ite)>= 0) ?  candidates[0].getSize() - coor[0]: coor[0];}while(ite++);
                                if (auto ite = candidates[0].getIterator()) do{(*ite) /=  candidates[0].getSize();}while(ite++);
                                exit(1);
                                break;
                            }
                        }else{

                            if (candidates[0].getSize() - torem.getSize() <= 1u){ // would leave nothing... try the highest second?
                                if (candidates[0].getSize() == torem.getSize()) {printf("Impossibleble!\n"); exit(1);}/*{ // wow, rounding error causes this
                                    acoor[0] = 0;
                                    for(acoor[1]=1;acoor[1] < torem.getSize();acoor[1]++) if (candidates[0][torem[acoor[0]]] < candidates[0][torem[acoor[1]]]) acoor[0] = acoor[1];
                                    torem.pop_swap(acoor[0]);
                                }*/

                                acoor[0] = 0;
                                for(acoor[1]=1;acoor[1] < torem.getSize();acoor[1]++) if (candidates[0][torem[acoor[0]]] < candidates[0][torem[acoor[1]]]) acoor[0] = acoor[1];
                                for(i=0;i<torem.getSize();i++) {
                                    if (i == acoor[0]) continue;
                                    j = candidates[0].find(torem[i]);
                                    if (j == 0xFFFFFFFF) {printf("erasing unexisting... !\n"); exit(1);}
                                    candidates[0].erase_from_iterator(j);
                                }
                                if (candidates[0].getSize() != 2) {printf("2 should remain. got %i\n", candidates[0].getSize()); exit(1);}
                                if (candidates[0].deref_key(0) == torem[acoor[0]]){
                                    candidates[0].deref(0) = 1.0;
                                    candidates[0].deref(1) = -1.0;
                                }else{
                                    candidates[0].deref(0) = -1.0;
                                    candidates[0].deref(1) = 1.0;
                                }
                                break;
                            }
                            for(i=0;i<torem.getSize();i++) {
                                j = candidates[0].find(torem[i]);
                                if (j == 0xFFFFFFFF) {printf("trying to remove un-existing!\n"); exit(1);}
                                candidates[0].erase_from_iterator(j);
                            }
                            sum[0] =0.0;
                            if (auto ite = candidates[0].getIterator()) do{sum[0] += *ite;}while(ite++);
                        }
                    }while(true);


                    if (auto ite = candidates[0].getIterator()){
                        while((*ite) > 0.0) ite++;
                        sum[0] = curcol[dir = ite()] / (*ite);

                        while(ite++){
                            if (*ite > 0.0) continue;
                            sum[1] = curcol[ite()] / (*ite);
                            if (sum[1] > sum[0]) {sum[0] = sum[1]; dir = ite();}
                        }
                        sum[0] = -sum[0]; // w has to be in 0-sum[0] range

                        tosolve_num2[0] = (target.rowhid_prior.mkInnerMult(candidates[0],curcol)+ target.rowhid_prior_constant) * ( 2.0 * ((tmptmp) ? trdata.data[coor[1]].getSize() : 1)  / prior_llterm[0]);
                        tosolve_num2[1] = d_scale[0];
                        tosolve_num2[2] = d_scale[1];
                        tosolve_den2.data[1] = 0.0; tosolve_den2.data[3] = 0.0;
                        tosolve_den2.data[2] = dd_scale.data[0];
                        tosolve_den2.data[5] = dd_scale.data[2];
                        tosolve_den2.data[4] = dd_scale.data[1];
                        tosolve_den2.data[0] = (target.rowhid_prior.mkInnerMult(candidates[0])+ target.rowhid_prior_constant) * ( 2.0 * ((tmptmp) ? trdata.data[coor[1]].getSize() : 1) / prior_llterm[0]);

                        for(i=0;i<trdata.data[coor[1]].getSize();i++){
                            coor[0]= trdata.data[coor[1]].deref_key(i);
                            if (!ExOp::isValid(target.par_colscale[coor[0]])) continue;
                            coeff = rever_buf[i].mkInnerProd(candidates[0]);
                            tosolve_num2[0] += coeff[0] * dd_buf[i][1] + coeff[1] * dd_buf[i][2];
                            tosolve_den2.data[0] += dd_buf[i][3] * coeff[0] * coeff[0] + dd_buf[i][4] * coeff[1] * coeff[1] + dd_buf[i][5] * coeff[0] * coeff[1];
                            tosolve_den2.data[1] += dd_buf[i][3] * coeff[0] + dd_buf[i][5] * coeff[1] * 0.5;
                            tosolve_den2.data[3] += dd_buf[i][4] * coeff[1] + dd_buf[i][5] * coeff[0] * 0.5;
                        }

                        tosolve_den2_tmp = tosolve_den2;
                        tosolve_num2_tmp = tosolve_num2;
                        tosolve_den2.toEigenTransform(funabsoft);
                        tosolve_num2 = tosolve_den2 * tosolve_num2;
                        tosolve_num2 *= hiddenalpha;

                        expexp = tosolve_num2[0] * tosolve_num2_tmp[0] + tosolve_num2[1] * tosolve_num2_tmp[1] + tosolve_num2[2] * tosolve_num2_tmp[2];
                        expexp2 = tosolve_den2_tmp.data[0] * tosolve_num2[0] * tosolve_num2[0] * 0.5;
                        expexp2 += tosolve_den2_tmp.data[2] * tosolve_num2[1] * tosolve_num2[1] * 0.5;
                        expexp2 += tosolve_den2_tmp.data[5] * tosolve_num2[2] * tosolve_num2[2] * 0.5;
                        expexp2 += tosolve_den2_tmp.data[1] * tosolve_num2[0] * tosolve_num2[1];
                        expexp2 += tosolve_den2_tmp.data[3] * tosolve_num2[0] * tosolve_num2[2];
                        expexp2 += tosolve_den2_tmp.data[4] * tosolve_num2[1] * tosolve_num2[2];

                        if (tosolve_num2[0] <= 0.0){ // went the other way... just update scale
                            if (tosolve_num2[0] == 0.0){
                                candidates[0] = curcol;
                            }else{
                                candidates[0] = curcol;
                                dd_scale_tmp = dd_scale; d_scale_tmp = d_scale;
                                d_scale = dd_scale.toEigenTransform(funabsoft) * d_scale;
                                d_scale *= hiddenalpha;
                            }
                        }else if (tosolve_num2[0] >= sum[0]){ // too far!, find solution on boundary...
                            d_scale[0] = tosolve_num2_tmp[1] + sum[0] * tosolve_den2.data[1];
                            d_scale[1] = tosolve_num2_tmp[2] + sum[0] * tosolve_den2.data[3];
                            d_scale = dd_scale.toEigenTransform(funabsoft) * d_scale;
                            d_scale *= hiddenalpha;
                            nb_funky++;
                            candidates[1] = candidates[0];
                            if (auto ite = candidates[0].getIterator()) do{
                                (*ite) = (*ite) * sum[0]  + (((j = curcol.find(ite())) == 0xFFFFFFFF) ? 0.0: curcol.deref(j));
                            }while(ite++);

                            if (auto ite = curcol.getIterator()) do{
                                if (candidates[0].find(ite())== 0xFFFFFFFF) candidates[0][ite()] = *ite;
                            }while(ite++);

                            j = candidates[0].find(dir);
                            if (j == 0xFFFFFFFF) {printf("trying to remove un-existing!\n");exit(1);}
                            candidates[0].erase_from_iterator(j);
                            sum[1] =0.0;
                            if (auto ite = candidates[0].getIterator()) {
                                do{sum[1] += *ite;}while(ite++);
                                if (fabs(sum[1] - 1.0) > 0.1) {printf("bound proj dir %i c... hidden state norm is %e! %e is fact\n", dir, sum[1], sum[0]);
                                candidates[0].show();
                                candidates[1].show();
                                curcol.show();
                                exit(1);}
                            }else {printf("got empty candidate!\n"); exit(1);}
                            candidates[1].toMemfree();
                        }else{ // all set!
                            d_scale[0] = tosolve_num2[1];
                            d_scale[1] = tosolve_num2[2];
                            candidates[1] = candidates[0];
                            if (auto ite = candidates[0].getIterator()) do{
                                (*ite) = (*ite) * tosolve_num2[0] + (((j = curcol.find(ite())) == 0xFFFFFFFF) ? 0.0: curcol.deref(j));
                            }while(ite++);
                            if (auto ite = curcol.getIterator()) do{
                                if (candidates[0].find(ite())== 0xFFFFFFFF) candidates[0][ite()] = *ite;
                            }while(ite++);
                            sum[1]= 0.0f;
                            if (auto ite = candidates[0].getIterator()) {
                                do{sum[1] += *ite;}while(ite++);
                                if (fabs(sum[1] - 1.0) > 0.1) {printf("normal proj rowr... hidden state norm is %e! %e is fact\n", sum[1], tosolve_num2[0]);
                                candidates[1].show();
                                candidates[0].show();
                                curcol.show();
                                LFH_ALIVE; exit(1);}
                            }else {printf("got empty candidate!\n"); exit(1);}
                            candidates[1].toMemfree();
                        }
                        candidates[0].makeMaximumAsFirst();
                        if ((candidates[0].deref_key(0) != curcol.deref_key(0))&&(!ignore_drop)){
                            // major change, add *last* boundary on the search space!
                            tosolve_num2[2] = ((i = curcol.find(candidates[0].deref_key(0))) == 0xFFFFFFFF) ? 0.0 : curcol.deref(i);
                            if (auto ite = curcol.getIterator()) {
                                tosolve_num2[0] = (*ite) - tosolve_num2[2]; // guarranied to be positive
                                tosolve_num2[1] = candidates[0].deref(0) -  (((i = candidates[0].find(ite())) == 0xFFFFFFFF) ? 0.0 : candidates[0].deref(i));
                                sum[1] = tosolve_num2[1] / (tosolve_num2[0] + tosolve_num2[1]);
                                j = 0;
                                while(ite++){
                                    if ((tosolve_num2[0] =  (*ite) - tosolve_num2[2]) <= 0) continue;
                                    tosolve_num2[1] = candidates[0].deref(0) -  (((i = candidates[0].find(ite())) == 0xFFFFFFFF) ? 0.0 : candidates[0].deref(i));
                                    sum[0] = tosolve_num2[1] / (tosolve_num2[0] + tosolve_num2[1]);
                                    if (sum[0] < sum[1]) {j = ite.i; sum[1] = sum[0];}
                                }
                                sum[0] = candidates[0].deref(0) * (1.0 - sum[1]) + sum[1] * (((i = curcol.find(candidates[0].deref_key(0))) == 0xFFFFFFFF) ? 0.0 : curcol.deref(i));
                                candidates[1][curcol.deref_key(j)] = sum[0];
                                candidates[1][candidates[0].deref_key(0)] = sum[0];
                                if (auto ite = curcol.getIterator()) do {
                                    if (candidates[1].find(ite()) == 0xFFFFFFFF) candidates[1][ite()] = (*ite) * sum[1] + (1.0 - sum[1]) * (((i = candidates[0].find(ite())) == 0xFFFFFFFF) ? 0.0 : candidates[0].deref(i));
                                }while(ite++);
                                if (auto ite = candidates[0].getIterator()) do {
                                    if (candidates[1].find(ite()) == 0xFFFFFFFF) candidates[1][ite()] = (*ite) * (1.0 - sum[1]);
                                }while(ite++);
                            }

                            if (auto ite = candidates[1].getIterator()) {
                                sum[0] =0.0;
                                do{
                                    sum[0] += *ite;
                                }while(ite++);
                                if (fabs(sum[0] - 1.0) > 0.1){printf("hidded state norm cancan1 is %e!\n", sum[0]);
                                    candidates[1].show();
                                    printf("frac %e\n", sum[1]);
                                    candidates[0].show();
                                    curcol.show();
                                LFH_ALIVE;exit(1);}
                            }else {printf("got empty candidate!\n"); exit(1);}

                        }
                    }else{
                        // stuck in corner again...
                        candidates[0] = curcol;
                        dd_scale_tmp = dd_scale; d_scale_tmp = d_scale;
                        d_scale = dd_scale.toEigenTransform(funabsoft) * d_scale;
                        d_scale *= hiddenalpha;
                    }
                }
                sum[0] = 0.0f;
                if (auto ite = candidates[0].getIterator()) {
                    do{
                        sum[0] += *ite;
                    }while(ite++);
                    if (fabs(sum[0] - 1.0) > 0.1){printf("hidded state norm is %e!\n", sum[0]);
                        candidates[0].show();
                        candidates[1].show();
                    LFH_ALIVE;exit(1);}
                }else {printf("got empty candidate!\n"); exit(1);}
                double tmp;

                //printf("cur coor %i/%i\n", coor[1], target.colhid.getNBcols());fflush(stdout);

                coeff = target.par_rowscale[coor[1]];
                coeff += d_scale;
                if ((tmp = (coeff[0]*coeff[0] + coeff[1]*coeff[1])) > 100.0){ // too far, send to boundary
                    tmp = pow(tmp, -0.5);
                    d_scale[0] = coeff[0] * tmp;
                    d_scale[1] = coeff[1] * tmp;
                    d_scale -= target.par_rowscale[coor[1]];
                }

                insist =0;
                while(true){

                    // Step 4 Evaluate Quantidates
                    for(j=0;j<2;j++){
                        future_statistics[j].toZero();future_statistic_trix[j].toZero();
                        curproj[j] = target.par_RM.mkBackMult(candidates[j]);
                        lls[j] = log(target.rowhid_prior.mkInnerMult(candidates[j])+ target.rowhid_prior_constant) * ((tmptmp) ? trdata.data[coor[1]].getSize() : 1) ;
                        if ((candidates[1].getSize() == 0)||(insist > 0)) break;
                    }
                                 //   printf("loopyb\n"); fflush(stdout);
                    for(i=0;i<trdata.data[coor[1]].getSize();i++){
                        coor[0]= trdata.data[coor[1]].deref_key(i);
                        if (!ExOp::isValid(target.par_colscale[coor[0]])) continue; // can do better...
                        for(j=0;j<2;j++){
                            coeff = curproj[j].mkInnerProd(target.colhid.getColumn(coor[0]));
                            coeff += target.par_colscale[coor[0]];
                            coeff += target.par_rowscale[coor[1]];
                            coeff += d_scale;
                            deflatedNBterm_wr_f_d2(ld.evilterms, trdata.data[coor[1]].deref(i), coeff[0], coeff[1]);
                           if (auto ite = target.colhid.getColumn(coor[0]).getIterator()) do{
                                future_statistics[j][ite()][0] += ld.evilterm_alt.d[0] * (*ite);
                                future_statistics[j][ite()][1] += ld.evilterm_alt.d[1] * (*ite);

                                trix_input[0] = ld.evilterms[3] * (*ite);
                                trix_input[1] = ld.evilterms[4] * (*ite);
                                trix_input[2] = ld.evilterms[5] * (*ite) * 0.5;
                                if (auto ite2 = target.colhid.getColumn(coor[0]).getIterator()) do{
                                    if (ite2() > ite()) continue;
                                    future_statistic_trix[j].data[ite2() + ((ite() * (ite() + 1)) >> 1)] += trix_input * (*ite2);
                                }while(ite2++);

                            }while(ite++);

                            lls[j] += ld.evilterms[0];
                            if ((candidates[1].getSize() == 0)||(insist > 0)) break;
                        }
                    }
                    //      printf("loopyb done\n"); fflush(stdout);

                    /*if (!checksafe) nb_funky++;
                    else {ld.evilterms[1]= (1.0 - (lls[0] - lls[3]) / prediction); accuracy += ld.evilterms[1]*ld.evilterms[1];}*/


                    // Step 4.1 Append LL from dropouts  TODO, ignore drop out for now
                    /*if (haschange){
                        for(i=0;i< nb_maxrowh.getSize();i++{


                        }
                    }*/

                    // Step 5 select best or do not update
                   // tmp = d_scale_tmp[0] * d_scale[0] + d_scale_tmp[1] * d_scale[1];
                   // tmp += (dd_scale_tmp.data[0] * d_scale[0] * d_scale[0] + dd_scale_tmp.data[2] * d_scale[1] * d_scale[1])* 0.5 + dd_scale_tmp.data[1] * d_scale[0] * d_scale[1];
                   // printf("%e incr (vs%e)\n", lls[0] - lls[3], tmp);

                    if (candidates[1].getSize() == 0){
                         j = (lls[0] > lls[2]) ? 0 : 2;
                    }else{
                        // worrying about dropouts!
                        lld.toZero();
                        acoor[0] = curcol.deref_key(0);
                        for(acoor[1] =0;acoor[1]< target.par_RM.sizes[1]; acoor[1]++){coor[0] = acoor[1];
                            lld[2] += row_nbmax_hcol(coor) * log(target.par_D(acoor));
                            lld[2] += (nb_maxcolh[coor[0]] - row_nbmax_hcol(coor)) * log(1.0 - target.par_D(acoor));
                        }
                        acoor[0] = candidates[0].deref_key(0);
                        for(acoor[1] =0;acoor[1]< target.par_RM.sizes[1]; acoor[1]++){coor[0] = acoor[1];
                            lld[0] += row_nbmax_hcol(coor) * log(target.par_D(acoor));
                            lld[0] += (nb_maxcolh[coor[0]] - row_nbmax_hcol(coor)) * log(1.0 - target.par_D(acoor));
                        }
                        if (candidates[1].deref_key(0) == curcol.deref_key(0)) lld[1] = lld[2];
                        else{
                            acoor[0] = candidates[1].deref_key(0);
                            for(acoor[1] =0;acoor[1]< target.par_RM.sizes[1]; acoor[1]++){coor[0] = acoor[1];
                                lld[1] += row_nbmax_hcol(coor) * log(target.par_D(acoor));
                                lld[1] += (nb_maxcolh[coor[0]] - row_nbmax_hcol(coor)) * log(1.0 - target.par_D(acoor));
                            }
                        }
                        j = (lls[0] + lld[0] > lls[2] + lld[2]) ? 0 : 2;
                        if ((lls[1] + ((lld[0] > lld[1]) ? lld[0]  : lld[1])) > lls[j] + lld[j]){
                            if (lld[0] > lld[1]) candidates[1].swapEntries(0,1);
                            j=1;
                        }
                        if ((j != 2)&&(candidates[j].deref_key(0) != curcol.deref_key(0))){
                            // registering change...
                            parascope[thrID].nbmax_diff[candidates[j].deref_key(0)]++;
                            parascope[thrID].nbmax_diff[curcol.deref_key(0)]--;
                            for(i = 0; i < trdata.data[coor[1]].getSize();i++){
                                acoor[1] = target.colhid.data[trdata.data[coor[1]].deref_key(i)].deref_key(0);
                                acoor[0] = curcol.deref_key(0);
                                parascope[thrID].nonzero_diff(acoor)--;
                                acoor[0] = candidates[j].deref_key(0);
                                parascope[thrID].nonzero_diff(acoor)++;
                            }
                        }
                    }

                    if (j != 2) break;
                    if (insist++ > 10) break;

                    if (auto ite = curcol.getIterator()) do{
                        candidates[0][ite()] += *ite;
                    }while(ite++);

                    if (auto ite = candidates[0].getIterator()) do{
                        *ite *= 0.5;
                    }while(ite++);
                    d_scale *= 0.5;


                    candidates[0].makeMaximumAsFirst();
                    if (candidates[0].deref_key(0) == curcol.deref_key(0)) candidates[1].toMemfree();

                }
                parascope[thrID].nb_insist += insist;
        //printf("curss coor %i/%i\n", coor[1], target.colhid.getNBcols());fflush(stdout);

                // Step 6 save scope for parameter update
                if (j == 2) { // failed to improve...
                    //printf("row %i failed %e %e\n", coor[1], target.par_rowscale[coor[1]][0], target.par_rowscale[coor[1]][1] );
                    candidates[0] = curcol;
                    newhidden.memmoveColumn(coor[1], candidates[0]);
                    newscale[coor[1]] = target.par_rowscale[coor[1]];
                    parascope[thrID].nb_fail++;
                }else{
                    newhidden.memmoveColumn(coor[1], candidates[j]);
                    newscale[coor[1]] = target.par_rowscale[coor[1]] + d_scale;
                }
                parascope[thrID].ll[1] += lls[j];
                candidates[1].toMemfree();

                future_factor.toMemfree();
                if (auto ite = newhidden.data[coor[1]].getIterator()) do{
                    future_factor[ite()][0] = *ite; future_factor[ite()][1] = *ite;
                }while(ite++);
                parascope[thrID].dll_dMR2.toAddBackOuterProd(future_statistics[j], future_factor);


                uint32_t dashift = target.par_RM.sizes[0] * target.par_RM.sizes[1];
                uint32_t k;
                for(i = 0,k=0; i< target.par_RM.sizes[1];i++){
                    for(uint32_t l=0;l<=i;l++,k++){
                        if (auto ite = newhidden.data[coor[1]].getIterator()) do{
                            uint32_t m = ite() + i * target.par_RM.sizes[0];
                            double* currow_A = parascope[thrID].dd_massive.data + (((m * (m +1)) >> 1) + l * target.par_RM.sizes[0]);
                            double* currow_B = parascope[thrID].dd_massive.data + ((((m + dashift) * (m + dashift +1)) >> 1) + l * target.par_RM.sizes[0]);
                            if (auto ite2 = newhidden.data[coor[1]].getIterator()) do{
                                if (ite2() > ite()) continue;
                                double tmp = (*ite) * (*ite2);
                                currow_A[ite2()] += future_statistic_trix[j].data[k][0] * tmp;
                                currow_B[ite2()] += future_statistic_trix[j].data[k][2] * tmp;
                                currow_B[ite2() + dashift] += future_statistic_trix[j].data[k][1] * tmp;
                            }while(ite2++);
                        }while(ite++);
                    }
                }
            }
        }
        lls[0] =0;
        for(coor[1]=colrange[thrID];coor[1] < colrange[thrID+1];coor[1]++){
            lls[0] += log(target.colhid_prior.mkInnerMult(target.colhid.getColumn(coor[1]))+ target.colhid_prior_constant) * ((tmptmp) ? data.data[coor[1]].getSize() : 1);
        }
        parascope[thrID].ll[0] += lls[0];
        parascope[thrID].ll[1] += lls[0];
    }
    if (parascope[thrID].use_stdout) {
        pbar.finish();
        printf("%e accuracy %c %e \n", parascope[thrID].alpha, (istimefor_columns) ? 'c': 'r', accuracy / ((istimefor_columns) ? colrange[thrID+1] - nb_funky : rowrange[thrID+1] - nb_funky));
    }
    return 1000;
}


};

 //   DataGrid<double, 2u> proj_X;
 //   Trianglix<double> row_w_cumul; row_w_cumul.setSize(par_RM.sizes[0]); // sum_i ZZ^t / Z^tUZ
 //   Trianglix<double> col_w_cumul; col_w_cumul.setSize(colprc.getSize()); // sum_j HH^t / H^tVH
   /* Tuple<double> meanproject;
    Tuple<double> rowsums[2];
    Tuple<double> colsums[2];
    TMatrix<double> Msums[2];
    Tuple<double> Msumstmp;
    Tuple<double> rowdenum[2];
    Tuple<double> coldenum[2];

    uint32_t cached_index;
    for(cached_index =0; cached_index <2u;cached_index++){
        rowsums[cached_index].setSize(data.getNBrows());
        colsums[cached_index].setSize(data.getNBcols());
        rowdenum[cached_index].setSize(data.getNBrows());
        coldenum[cached_index].setSize(data.getNBcols());
        Msums[cached_index].setSizes(par_RM.sizes[0],par_RM.sizes[1]);
    }cached_index =0;
    Msumstmp.setSize(par_RM.sizes[0]);*/

    Task scp(*this,data,excl_list);

    scp.ignore_drop = false;
    scp.hiddenalpha = 1.0;

    double that_ll_constant=0.0;
    if (auto ite = data.getIterator()) do{
        that_ll_constant -= lngamma(1.0 + *ite);
    }while(ite++);

    Tuple<uint32_t, 2u> coor, coor2, coorC, coorR;
    coor[0] = par_RM.sizes[0];
    coor[1] = par_RM.sizes[1];
    TMatrix< Tuple<double ,2u> > oldRM;
    double old_LL = -DBL_MAX;

    uint32_t i,j,k;
    scp.parascope.setSize(nbthread);
    scp.rowrange.setSize(nbthread+1);
    scp.colrange.setSize(nbthread+1);
    scp.rowrange[0] =0; scp.colrange[0] =0;
    for(i=0;i<nbthread;i++) {
        scp.parascope[i].setSizes(coor,0);
        //parascope[i].setPointers((Tuple<double>*)rowdenum,(Tuple<double>*)coldenum, &cached_index);
        scp.rowrange[i+1] = ((i+1) * scp.trdata.getNBcols()) / nbthread;
        scp.colrange[i+1] = ((i+1) * data.getNBcols()) / nbthread;
        scp.parascope[i].alpha = 0.001;
        scp.parascope[i].use_stdout = (i == nbthread-1);
    }

    ProgressBarPrint pbar(20);
    double tmp,tmpsum, ll,ll_best;
    FunctionScope<double&,double&,double> funabsoft(ExCo<double>::toAbsoftInverse,pow(0.5, 256.0f));
    Trianglix<double,2u> single_denum;
    Tuple<double, 2u> single_num;
    Tuple<double,0u,TUPLE_FLAG_REMOTE_MEMORY> dainput;

    uint32_t rowite;
    int step;
    CurvatureSearchScope cs;
    //scp.istimefor_columns = !scp.istimefor_columns;
    uint32_t fuz =0;
    double expected,expected2,dropLL,last_dll,rectsize;
    Tuple<double,0u> toproj; toproj.setSize(2 * par_RM.sizes[0] * par_RM.sizes[1]);
    Tuple<double,0u> ddproj; ddproj.setSize(2 * par_RM.sizes[0] * par_RM.sizes[1]);
    cs.init(2 * par_RM.sizes[0] * par_RM.sizes[1], 1.0);
    double logalpha = -1.0;
    dropLL = 0.0f;
    if (auto ite = par_D.getIterator()) do{
        (*ite) = ((double)scp.maxgrid_nonzero_count(ite())) / (scp.nb_maxrowh[ite()[0]] * scp.nb_maxcolh[ite()[1]]); // p = c_ij / (s_i*t_j)
    }while(ite++);

    while(true){
        scp.newhidden.setNBcols((scp.istimefor_columns)? colhid.getNBcols() : rowhid.getNBcols());
        scp.newscale.setSize(scp.newhidden.getNBcols());
 //       scp.new_nbmax_h.setSizes((scp.istimefor_columns)? scp.row_nbmax_hcol.getDims() : scp.col_nbmax_hrow.getDims()  );

        scp.tmptmp = false; //(fuz >= (nbstep/2)); // heavy/light hidden prior
        //scp.heuristictime = true;

        if ((fuz % 10) == 5) doRecenter(); // hehe

        for(i=nbthread-1;i>0;i--) tb.submit(scp,i);
        tb.submit_ThenWait(scp, 0);
        for(i=1;i<nbthread;i++) {
            ExOp::toAdd(scp.parascope[0].ll,scp.parascope[i].ll);
            scp.parascope[0].alpha += scp.parascope[i].alpha;
            scp.parascope[0].nonzero_diff += scp.parascope[i].nonzero_diff;
            scp.parascope[0].nbmax_diff += scp.parascope[i].nbmax_diff;
            scp.parascope[0].nb_fail += scp.parascope[i].nb_fail;
            scp.parascope[0].nb_insist += scp.parascope[i].nb_insist;
        }
        if (!scp.ignore_drop) {
            scp.parascope[0].nonzero_diff += scp.maxgrid_nonzero_count;
            dropLL =0.0;
            if (scp.istimefor_columns){
                for(i=0;i<par_RM.sizes[1];i++) scp.parascope[0].nbmax_diff[i] += scp.nb_maxcolh[i];
                if (auto ite = par_D.getIterator()) do{
                    rectsize = ((scp.nb_maxrowh[ite()[0]] * scp.parascope[0].nbmax_diff[ite()[1]]));
                    if ((scp.parascope[0].nonzero_diff(ite()) != 0)&&(scp.parascope[0].nonzero_diff(ite()) != scp.nb_maxrowh[ite()[0]] * scp.parascope[0].nbmax_diff[ite()[1]])){
                        dropLL += log((double)scp.parascope[0].nonzero_diff(ite())) * scp.parascope[0].nonzero_diff(ite());
                        dropLL += log((double)(rectsize - scp.parascope[0].nonzero_diff(ite()))) * (rectsize - scp.parascope[0].nonzero_diff(ite()));
                        dropLL -= log((double)(rectsize)) * rectsize;
                        dropLL += lngamma((double)rectsize);
                        dropLL -= lngamma((double)(rectsize- scp.parascope[0].nonzero_diff(ite())));
                        dropLL -= lngamma((double)(scp.parascope[0].nonzero_diff(ite())));
                    }// else iszero
                }while(ite++);

            }else{
                for(i=0;i<par_RM.sizes[0];i++) scp.parascope[0].nbmax_diff[i] += scp.nb_maxrowh[i];
                if (auto ite = par_D.getIterator()) do{
                    rectsize = ((double)(scp.parascope[0].nbmax_diff[ite()[0]] * scp.nb_maxcolh[ite()[1]]));
                    if ((scp.parascope[0].nonzero_diff(ite()) != 0)&&(scp.parascope[0].nonzero_diff(ite()) != scp.parascope[0].nbmax_diff[ite()[0]] * scp.nb_maxcolh[ite()[1]])){
                        dropLL += log((double)scp.parascope[0].nonzero_diff(ite())) * scp.parascope[0].nonzero_diff(ite());
                        dropLL += log((double)(rectsize - scp.parascope[0].nonzero_diff(ite()))) * (rectsize - scp.parascope[0].nonzero_diff(ite()));
                        dropLL -= log(rectsize) * rectsize;
                        dropLL += lngamma(rectsize);
                        dropLL -= lngamma((double)(rectsize- scp.parascope[0].nonzero_diff(ite())));
                        dropLL -= lngamma((double)(scp.parascope[0].nonzero_diff(ite())));
                    }// else iszero
                }while(ite++);
            }
            //printf("dropoutLL %e -> %e\n", last_dll, dropLL);
            last_dll = dropLL;
            scp.parascope[0].ll[0] += dropLL; // not exact (dropout), but not important
            scp.parascope[0].ll[1] += dropLL;
        }
        scp.parascope[0].alpha  /= nbthread;
        for(i=1;i<nbthread;i++) scp.parascope[i].alpha = scp.parascope[0].alpha ;

        rectsize = ((double)scp.parascope[0].nb_insist) / ((scp.istimefor_columns) ? scp.data.getNBcols() : scp.trdata.getNBcols());
        if (rectsize > 0.5) scp.hiddenalpha *= exp(-1);
        else scp.hiddenalpha = exp( 0.75 * log(scp.hiddenalpha) -0.25 * rectsize);
        rectsize = ((scp.istimefor_columns) ? scp.data.getNBcols() : scp.trdata.getNBcols());

        printf("%i Likelihood %e -> %e (%c In rate%e FF %e)\n", fuz, scp.parascope[0].ll[0] + that_ll_constant, scp.parascope[0].ll[1] + that_ll_constant, scp.istimefor_columns ? 'C' : 'R', ((double)scp.parascope[0].nb_insist) / rectsize, ((double)scp.parascope[0].nb_fail) / rectsize);
        printf("expected change %e, got %e  (ddmag %e)\n", expected + expected2, scp.parascope[0].ll[0] - old_LL, expected2);

        if ((old_LL < scp.parascope[0].ll[1])||(fuz == 0)||((scp.tmptmp)&&(fuz == (nbstep/2)))){
            expected = 0.0;
            expected2 = 0.0;
            oldRM = par_RM;
            old_LL = scp.parascope[0].ll[1];
            if (scp.istimefor_columns){
                for(i=0;i< colhid.getNBcols();i++){
                    if ((colhid.data[i].getSize() == 0)||(scp.newhidden.data[i].getSize() == 0)) continue;
                    if (colhid.data[i].deref_key(0) != scp.newhidden.data[i].deref_key(0)){
                        for(j = 0; j < data.data[i].getSize();j++){
                            coor[1] = data.data[i].deref_key(j);
                            coor[0] = colhid.data[i].deref_key(0);
                            scp.row_nbmax_hcol(coor)--;
                            coor[0] = scp.newhidden.data[i].deref_key(0);
                            scp.row_nbmax_hcol(coor)++;
                        }
                    }
                }
                for(i=0;i<par_RM.sizes[1];i++) scp.nb_maxcolh[i] = scp.parascope[0].nbmax_diff[i];
                colhid.toMemmove(scp.newhidden);
                par_colscale.toMemmove(scp.newscale);
             }else{
                for(i=0;i< rowhid.getNBcols();i++){
                    if ((rowhid.data[i].getSize() == 0)||(scp.newhidden.data[i].getSize() == 0)) continue;
                    if (rowhid.data[i].deref_key(0) != scp.newhidden.data[i].deref_key(0)){
                        for(j = 0; j < scp.trdata.data[i].getSize();j++){
                            coor[1] = scp.trdata.data[i].deref_key(j);
                            coor[0] = rowhid.data[i].deref_key(0);
                            scp.col_nbmax_hrow(coor)--;
                            coor[0] = scp.newhidden.data[i].deref_key(0);
                            scp.col_nbmax_hrow(coor)++;
                        }
                    }
                }
                for(i=0;i<par_RM.sizes[0];i++) scp.nb_maxrowh[i] = scp.parascope[0].nbmax_diff[i];
                rowhid.toMemmove(scp.newhidden);
                par_rowscale.toMemmove(scp.newscale);
            }

            for(i=1;i<nbthread;i++) {
                scp.parascope[0].dd_massive += scp.parascope[i].dd_massive;
                scp.parascope[0].dll_dMR2 += scp.parascope[i].dll_dMR2;
            }
            // uptate D, close close form hurray!
            scp.maxgrid_nonzero_count = scp.parascope[0].nonzero_diff;
            if (auto ite = par_D.getIterator()) do{
                (*ite) = ((double)scp.maxgrid_nonzero_count(ite())) / (scp.nb_maxrowh[ite()[0]] * scp.nb_maxcolh[ite()[1]]);
            }while(ite++);
            if (fuz >= nbstep) break;

            // fill symmetries... yea
            if (auto ite = scp.parascope[0].dd_massive.getIterator()) do{
                if ((ite()[0] % par_RM.sizes[0]) > (ite()[1] % par_RM.sizes[0])){
                    (*ite) = scp.parascope[0].dd_massive.cell( ((ite()[1] % par_RM.sizes[0]) + ((ite()[0] / par_RM.sizes[0]) * par_RM.sizes[0])) , ((ite()[0] % par_RM.sizes[0]) + ((ite()[1] / par_RM.sizes[0]) * par_RM.sizes[0])) );
                }
            }while(ite++);
            if (auto ite = scp.parascope[0].dd_massive.getIterator()) do{
                if ((ite()[0] % (par_RM.sizes[0] * par_RM.sizes[1])) > (ite()[1] % (par_RM.sizes[0] * par_RM.sizes[1]))){
                    (*ite) = scp.parascope[0].dd_massive.cell(ite()[1] % (par_RM.sizes[0] * par_RM.sizes[1]), (par_RM.sizes[0] * par_RM.sizes[1]) + ite()[0] % (par_RM.sizes[0] * par_RM.sizes[1]));
                }
            }while(ite++);

            i=0;
            if (auto ite = scp.parascope[0].dll_dMR2.getIterator()) do{
                toproj[i + par_RM.sizes[0] * par_RM.sizes[1]] = (*ite)[1];
                toproj[i++] = (*ite)[0];
            }while(ite++);
            scp.parascope[0].dd_massive.toEigenTransform(funabsoft);
            toproj = scp.parascope[0].dd_massive * toproj;


            #ifdef Rcpp_hpp
          //  if (scp.istimefor_columns){
            i=0;
                if (auto ite = par_RM.getIterator()) do{
                    expected += exp(logalpha) * toproj[i] * scp.parascope[0].dll_dMR2(ite())[0];
                    expected += exp(logalpha) * toproj[i + par_RM.sizes[0] * par_RM.sizes[1]] * scp.parascope[0].dll_dMR2(ite())[1];
                    (*ite)[1] += exp(logalpha) * toproj[i + par_RM.sizes[0] * par_RM.sizes[1]];
                    (*ite)[0] += exp(logalpha) * toproj[i++];
                }while(ite++);
           // }

            #endif

            scp.istimefor_columns = !scp.istimefor_columns;
            logalpha = (logalpha * 127.0) / 128.0;

        }else{ // failed... fall backward and do not update
            if (fuz >= nbstep + 10) break;
            printf("was tried:\n");
            par_RM.show();
            printf("backup: %e\n", logalpha);
            oldRM.show();
            printf("FAILFAILFAILFAILFAILFAILFAILFAILFAILFAILFAIL %e\n", old_LL + that_ll_constant); fflush(stdout);
            par_RM = (oldRM * 3.0 + par_RM) * 0.25;
            printf("FAILFAILFAILFAILFAILFAILFAILFAILFAILFAILFAIL %e\n", scp.parascope[0].ll[0] + that_ll_constant); fflush(stdout);
            expected *= 0.25;
            logalpha -= 1.0;
        }
        fuz++;
     }



     if (scp.ignore_drop){ // update this, was not used
         scp.maxgrid_nonzero_count.toZero();
         if (auto ite = data.getIterator()) do{
            coor[0] = rowhid.data[ite()[0]].deref_key(0);
            coor[1] = colhid.data[ite()[1]].deref_key(0);
            scp.maxgrid_nonzero_count(coor)++;
         }while (ite++);
         scp.nb_maxrowh.toZero();
         scp.nb_maxcolh.toZero();
         for(coor[1]=0;coor[1]<rowhid.getNBcols();coor[1]++) {
            rowhid.makeMaximumAsFirst(coor[1]);
            scp.nb_maxrowh[rowhid.data[coor[1]].deref_key(0)]++;
        }
        for(coor[1]=0;coor[1]<colhid.getNBcols();coor[1]++){
            colhid.makeMaximumAsFirst(coor[1]);
            scp.nb_maxcolh[colhid.data[coor[1]].deref_key(0)]++;
        }
         if (auto ite = par_D.getIterator()) do{
            (*ite) = ((double)scp.maxgrid_nonzero_count(ite())) / (scp.nb_maxrowh[ite()[0]] * scp.nb_maxcolh[ite()[1]]);
         }while(ite++);
    }


    doRecenter();
}

LFHTEMP ERRCODE BiClassifier<C,R,S>::run2D_EM_dropout(ThreadBase &tb, Tuple< double > &fout_rowdrop, Tuple< double > &fout_coldrop, const SparseMatrix<C>& data, uint32_t nbstep){
class Task : public Event<uint32_t>{
public:
ThreadBase &tb;
const SparseMatrix<C>& input;
SparseMatrix<C> trinput;
Tuple<double> &fout_rowdrop;
Tuple<double> &fout_coldrop;
Tuple<double> row_nextguess;
Tuple<double> col_nextguess;
Vector< KeyElem<double, uint32_t> > freq_rows;
Vector< KeyElem<double, uint32_t> > freq_cols;
// Tuple<Tuple<double ,3u> > sufstats;
Tuple<uint32_t> row_range;
Tuple<uint32_t> col_range;
myHashmap<uint32_t, void> finite_row;
myHashmap<uint32_t, void> finite_col;
uint32_t* aprogr;
uint32_t nbnonzero;
double LL_buffer;
double LL_buffer_chk;
uint32_t nbstep;
uint32_t tasks;
uint32_t thrprint;
Task(ThreadBase &_tb, const SparseMatrix<C>&  _input, Tuple<double>& _rowdrop, Tuple<double>& _coldrop,uint32_t _nbstep):
    tb(_tb),input(_input), fout_rowdrop(_rowdrop), fout_coldrop(_coldrop),nbstep(_nbstep){}
operator bool(){return (*this)() != 0;}
operator ERRCODE (){return (*this)();}
ERRCODE operator()(){
    fout_rowdrop.setSize(input.computeNBrows());
    fout_coldrop.setSize(input.getNBcols());

    uint32_t nbthread = tb.getSize();
    uint32_t i,j;

    row_range.setSize(nbthread+1);
    col_range.setSize(nbthread+1);

    row_nextguess.setSize(fout_rowdrop.getSize());
    col_nextguess.setSize(fout_coldrop.getSize());

    // sufstats.setSize( (fout_rowdrop.getSize() > fout_coldrop.getSize()) ? fout_rowdrop.getSize(): fout_coldrop.getSize());

    Vector< KeyElem<uint32_t, uint32_t> > freq_rows, freq_cols;
    trinput = input.mkTranspose();
    freq_rows.setSize(fout_rowdrop.getSize());
    for(i=0;i<fout_rowdrop.getSize();i++){
        freq_rows[i].d = i;
        freq_rows[i].k = trinput.data[i].getSize();
    }
    freq_rows.sort();
    freq_cols.setSize(fout_coldrop.getSize());
    for(i=0;i<fout_coldrop.getSize();i++){
        freq_cols[i].d = i;
        freq_cols[i].k = input.data[i].getSize();
    }
    freq_cols.sort();

    uint32_t totflt, tmpflt;
    uint32_t rang[2][2];
    rang[0][0] = 0; // nb emptyrow
    rang[1][0] = 0; // nb emptycol
    rang[0][1] = 1; // nb fillrow + 1
    rang[1][1] = 1; // nb fillcol + 1
    tmpflt =0;
    bool breatime= false;
    do{
        totflt = tmpflt;
        while((freq_rows[rang[0][0]].k) < rang[1][1]) {
            rang[0][0]++;
            if ((rang[0][0] + rang[0][1]) > fout_rowdrop.getSize()) {breatime = true; break;}
        }
        if (breatime) break;
        while((freq_rows[fout_rowdrop.getSize() - rang[0][1]].k) >= fout_coldrop.getSize() - rang[1][0]) {
            rang[0][1]++;
            if ((rang[0][0] + rang[0][1]) > fout_rowdrop.getSize()) {breatime = true; break;}
        }
        if (breatime) break;
        while((freq_cols[rang[1][0]].k) < rang[0][1]) {
            rang[1][0]++;
            if ((rang[1][0] + rang[1][1]) > fout_coldrop.getSize()) {breatime = true; break;}
        }
        if (breatime) break;
        while((freq_cols[fout_coldrop.getSize() - rang[1][1]].k) >= fout_rowdrop.getSize() - rang[0][0]) {
            rang[1][1]++;
            if ((rang[1][0] + rang[1][1]) > fout_coldrop.getSize())  {breatime = true; break;}
        }
        if (breatime) break;
        if ((rang[1][0] + rang[1][1]) > fout_coldrop.getSize())  break;
        tmpflt = rang[0][0] + rang[0][1] + rang[1][0] + rang[1][1];
    }while(tmpflt != totflt);

    printf("recursively filtered %i+%i rows and %i+%i columns (empty + filled)\n", rang[0][0], rang[0][1] -1, rang[1][0], rang[1][1] -1);

    uint32_t upper[2];
    upper[0] = fout_rowdrop.getSize() - rang[0][1] + 1;
    upper[1] = fout_coldrop.getSize() - rang[1][1] + 1;

    for(i = 0;i < rang[0][0];i++) ExCo<double>::toNegInfinity(fout_rowdrop[freq_rows[i].d]);
    for(i = upper[0];i < fout_rowdrop.getSize();i++) ExCo<double>::toInfinity(fout_rowdrop[freq_rows[i].d]);
    for(i = 0;i < rang[0][1];i++) ExCo<double>::toNegInfinity(fout_coldrop[freq_cols[i].d]);
    for(i = upper[1];i < fout_coldrop.getSize();i++) ExCo<double>::toInfinity(fout_coldrop[freq_cols[i].d]);

    if (((rang[1][0] + rang[1][1]) > fout_coldrop.getSize())||((rang[0][0] + rang[0][1]) > fout_rowdrop.getSize())) {
        // everything got filtered... how fun?
        printf("Everything got filtered!");
        for(i = rang[0][0]; i < upper[0];i++) ExCo<double>::toUndefined(fout_rowdrop[freq_rows[i].d]);
        for(i = rang[0][1]; i < upper[1];i++) ExCo<double>::toUndefined(fout_coldrop[freq_cols[i].d]);
        return 1;
    }
    j = (fout_coldrop.getSize() + 1 - rang[1][0] - rang[1][1]);

    for(i = rang[0][0]; i < upper[0];i++) {
        finite_row.addEntry(freq_rows[i].d);
        //printf("%i <- %i - %i / %i\n", freq_rows[i].d, freq_rows[i].k, rang[1][1], j);
        fout_rowdrop[freq_rows[i].d] = logit(((double)freq_rows[i].k - rang[1][1] + 1) / j);
        if (fabs(fout_rowdrop[freq_rows[i].d]) < 1.0) row_nextguess[freq_rows[i].d] = -fout_rowdrop[freq_rows[i].d];
        else row_nextguess[freq_rows[i].d] = fout_rowdrop[freq_rows[i].d] * 0.5;
    }

    j = (fout_rowdrop.getSize() + 1 - rang[0][0] - rang[0][1]);
    for(i = rang[1][0]; i < upper[1];i++) {
        finite_col.addEntry(freq_cols[i].d);
        //printf("%i <- %i - %i / %i\n", freq_cols[i].d, freq_cols[i].k, rang[0][1], j);
        fout_coldrop[freq_cols[i].d] = logit(((double)freq_cols[i].k - rang[0][1] + 1) / j);
        if (fabs(fout_coldrop[freq_cols[i].d]) < 1.0) col_nextguess[freq_cols[i].d] = -fout_coldrop[freq_cols[i].d];
        else col_nextguess[freq_cols[i].d] = fout_coldrop[freq_cols[i].d] * 0.5;
    }
    for(i=0;i<=nbthread;i++) {
        row_range[i] = (i * finite_row.getSize()) / nbthread;
        col_range[i] = (i * finite_col.getSize()) / nbthread;
    }

    // remove empty and filled rows / columns recursively!
    /*
    tasks = 2;aprogr =0; aprogr_factor = 1.0 / row_range[nbthread];
    thrprint =0; progb.start("Init dropout parameter guess for rows");
    for(i=nbthread-1;i>0;i--) tb.startEvent(this,i);
    tb.startEvent_ThenWait(this,i);
    progb.finish();
    tasks = 3;aprogr =0; aprogr_factor = 1.0 / col_range[nbthread];
    thrprint =0; progb.start("Init dropout parameter guess for cols");
    for(i=nbthread-1;i>0;i--) tb.startEvent(this,i);
    tb.startEvent_ThenWait(this,i);
    progb.finish();*/
    thrprint =0;
    for(j=0;j<10;j++){
        LL_buffer_chk =0;LL_buffer = 0.0; tasks = 0;tb.startProgress("Update Row Parameters", row_range[nbthread]);
        for(i=nbthread-1;i>0;i--) tb.startEvent(this,i);
        tb.startEvent_ThenWait(this,i);
        printf("New LL %e (old %e)\n", LL_buffer, LL_buffer_chk);
        LL_buffer_chk =0;LL_buffer = 0.0; tasks = 1;tb.startProgress("Update Col Parameters", col_range[nbthread]);
        for(i=nbthread-1;i>0;i--) tb.startEvent(this,i);
        tb.startEvent_ThenWait(this,i);
        printf("New LL %e (old %e)\n", LL_buffer, LL_buffer_chk);
    }
}
uint32_t operator()(uint32_t threadID){
    uint32_t itesel, sel, ort;
    FuzzyLineSearch fls;
    uint32_t k;
    switch(tasks){
    case 0:
        for(itesel = row_range[threadID]; itesel < row_range[threadID+1];itesel++){
            sel = finite_row.deref(itesel);

            fls.initGuesses(fout_rowdrop[sel],row_nextguess[sel]);
//                printf("c%e %e %e is guess\n", fls.guess[0],fls.guess[1],fls.guess[2]);
            for(ort = 0; ort < finite_col.getSize();ort++){
              //  printf("%i, %i: %e\n", sel, finite_col.deref(ort), fls.guess[0] + fout_coldrop[finite_col.deref(ort)] );
                if (input.data[finite_col.deref(ort)].find(sel) == 0xFFFFFFFF){
                    for(k=0;k<3;k++) fls.value[k] -= log_e_x_p1(fls.guess[k] + fout_coldrop[finite_col.deref(ort)]);
                }else{
                    for(k=0;k<3;k++) fls.value[k] -= log_e_x_p1(-fls.guess[k] - fout_coldrop[finite_col.deref(ort)]);
                }
            }
            // pick smallest of points, find best next guess
            LL_buffer += fls.solveArgmax(fout_rowdrop[sel],row_nextguess[sel]);
            LL_buffer_chk += fls.value[0];
            tb.updateProgress(threadID);
        }
        break;
    case 1:
        for(itesel = col_range[threadID]; itesel < col_range[threadID+1];itesel++){
            sel = finite_col.deref(itesel);
            fls.initGuesses(fout_coldrop[sel],col_nextguess[sel]);
//                printf("c%e %e %e is guess\n", fls.guess[0],fls.guess[1],fls.guess[2]);
            for(ort = 0; ort < finite_row.getSize();ort++){
//                    printf("c%e is a\n", fout_rowdrop[finite_row.deref(ort)]);
             //   printf("%i, %i: %e\n", finite_row.deref(ort), sel, fls.guess[0] + fout_rowdrop[finite_row.deref(ort)] );
                if (trinput.data[finite_row.deref_key(ort)].find(sel) == 0xFFFFFFFF){
                    for(k=0;k<3;k++) fls.value[k] -= log_e_x_p1(fls.guess[k] + fout_rowdrop[finite_row.deref(ort)]);
                }else{
                    for(k=0;k<3;k++) fls.value[k] -= log_e_x_p1(-fls.guess[k] - fout_rowdrop[finite_row.deref(ort)]);
                }
            }
            LL_buffer += fls.solveArgmax(fout_coldrop[sel],col_nextguess[sel]);
            LL_buffer_chk += fls.value[0];
            tb.updateProgress(threadID);
        }
        break;
   /* case 2: {
        for(sel = row_range[threadID]; sel < row_range[threadID+1];sel++, aprogr++){
            if (thrprint == threadID) progb.update(aprogr_factor * aprogr);
            fout_rowdrop[sel] = logit(((double)trinput.data.getSize()) / input.getNBcols());
            if (ExCo<double>::isInfinite(fout_rowdrop[sel])) infinite_row.addEntry(sel);
            else{
                if (fabs(fout_rowdrop[sel]) < 1.0) row_nextguess[sel] = -fout_rowdrop[sel];
                else row_nextguess[sel] = fout_rowdrop[sel] * 0.5;
            }
        }
        for(sel = col_range[threadID]; sel < col_range[threadID+1];sel++, aprogr++){
            if (thrprint == threadID) progb.update(aprogr_factor * aprogr);
            fout_coldrop[sel] = logit(((double)input.data.getSize()) / trinput.getNBcols());
            if (ExCo<double>::isInfinite(fout_coldrop[sel])) infinite_col.addEntry(sel);
            else{
                if (fabs(fout_coldrop[sel]) < 1.0) row_nextguess[sel] = -fout_coldrop[sel];
                else row_nextguess[sel] = fout_coldrop[sel] * 0.5;
            }
        }
    }break;*/
    }
return 0;}
}; // end of Task Scope
return (ERRCODE) Task(tb,data,fout_rowdrop,fout_coldrop,nbstep);}

LFHTEMP template<class P> void BiClassifier<C,R,S>::run2D_scaleEM(const SparseMatrix<C>& data, P& tb, uint32_t nbstep, unsigned int nbthread){

}
LFHTEMP void BiClassifier<C,R,S>::doRecenter(){
    Tuple<double > hidcol_frequency, hidrow_frequency;
    Tuple<Tuple<double,2> > row_average, col_average;
    Tuple<double,2> factaverage[3];
    Tuple<uint32_t,2> counts;
    hidrow_frequency.setSize(par_RM.sizes[0]).toZero();
    hidcol_frequency.setSize(par_RM.sizes[1]).toZero();
    row_average.setSize(par_RM.sizes[0]).toZero();
    col_average.setSize(par_RM.sizes[1]).toZero();
    ExOp::toZero(factaverage); counts.toZero();
    uint32_t i;
    for(i=0;i< par_rowscale.getSize();i++){
        if (!ExOp::isValid(par_rowscale[i][0])) continue;
        if (auto ite = rowhid.data[i].getIterator()) do{
            hidrow_frequency[ite()] += *ite;
        }while(ite++);
        counts[0]++;
    }
    hidrow_frequency /= counts[0];
    for(i=0;i< par_colscale.getSize();i++){
        if (!ExOp::isValid(par_colscale[i][0])) continue;
        if (auto ite = colhid.data[i].getIterator()) do{
            hidcol_frequency[ite()] += *ite;
        }while(ite++);
        counts[1]++;
    }
    hidcol_frequency /= counts[1];
    if (auto ite = par_RM.getIterator()) do{
        row_average[ite()[0]] += (*ite) * hidcol_frequency[ite()[1]];
        col_average[ite()[1]] += (*ite) * hidrow_frequency[ite()[0]];
    }while(ite++);

    for(i=0;i<par_RM.sizes[0];i++) factaverage[2] += row_average[i] * hidrow_frequency[i];
    for(i=0;i<par_RM.sizes[1];i++) factaverage[2] += col_average[i] * hidcol_frequency[i];
    factaverage[2] *= -0.25;
    for(i=0;i<par_RM.sizes[0];i++) row_average[i] += factaverage[2];
    for(i=0;i<par_RM.sizes[1];i++) col_average[i] += factaverage[2];

    if (auto ite = par_RM.getIterator()) do{
        (*ite) -= row_average[ite()[0]] + col_average[ite()[1]];
    }while(ite++);

    for(i=0;i< par_rowscale.getSize();i++){
        if (!ExOp::isValid(par_rowscale[i][0])) continue;
        if (auto ite = rowhid.data[i].getIterator()) do{
            par_rowscale[i] += row_average[ite()] * (*ite);
        }while(ite++);
        factaverage[0] += par_rowscale[i];
    }
    for(i=0;i< par_colscale.getSize();i++){
        if (!ExOp::isValid(par_colscale[i][0])) continue;
        if (auto ite = colhid.data[i].getIterator()) do{
            par_colscale[i] += col_average[ite()] * (*ite);
        }while(ite++);
        factaverage[1] += par_colscale[i];
    }
    factaverage[0] /= counts[0];
    factaverage[1] /= counts[1];
    factaverage[0] -= factaverage[1];
    factaverage[0] *= 0.5;
    for(i=0;i< par_rowscale.getSize();i++) if (ExOp::isValid(par_rowscale[i][0])) par_rowscale[i] -= factaverage[0];
    for(i=0;i< par_colscale.getSize();i++) if (ExOp::isValid(par_colscale[i][0])) par_colscale[i] += factaverage[0];

    // checkcheck
    row_average.toZero();
    col_average.toZero();
    if (auto ite = par_RM.getIterator()) do{
        row_average[ite()[0]] += (*ite) * hidcol_frequency[ite()[1]];
        col_average[ite()[1]] += (*ite) * hidrow_frequency[ite()[0]];
    }while(ite++);
    printf("did recenter!\n");
    row_average.show();
    col_average.show();
}

// unsigned char scaleconv(double d){return (unsigned char) (d * 256.0);}
LFHTEMP void BiClassifier<C,R,S>::exportTiff(const char* folderpath)const{
    Tuple<uint32_t,2u> coor;
    char buffer[1024];
    uint32_t i = strlen(folderpath);
    memcpy(buffer, folderpath, i);
    strcpy(buffer+i, "Hidden_rows.tif");
    TiffFile tfor(buffer,true);
    strcpy(buffer+i, "Hidden_cols.tif");
    TiffFile tfoc(buffer,true);
    strcpy(buffer+i, "Expectation.tif");
    TiffFile tfom(buffer,true);
    //strcpy(buffer+i, "Mat_var.tif");
    //TiffFile tfov(buffer,true);
//    ScaleConvert<unsigned char, double, double> sc(256.0);

    DataGrid<unsigned char, 2u> daout; // = ExOp::apply(sc,rowhid());
    DataGrid<float, 2u> daouf; // = ExOp::apply(sc,rowhid());

    coor[0] = par_RM.sizes[0];
    coor[1] = rowhid.getNBcols();
    daout.setSizes(coor);
    coor[0] = 2;
    daouf.setSizes(coor);
    unsigned char* optr = daout.data;
    double sum;
    for(coor[1]=0;coor[1] < daout.dims[1];coor[1]++){
        if (rowhid.data[coor[1]].getSize() == 0) printf("error! row %i has no state!\n", coor[1]);
        sum =0;
        for(coor[0]=0;coor[0] < daout.dims[0];coor[0]++){
            if ((i = rowhid.data[coor[1]].find(coor[0])) == 0xFFFFFFFF) *(optr++) = 0;
            else {
                *(optr++) = (unsigned char)(rowhid.data[coor[1]].deref(i) * 255);
                sum += rowhid.data[coor[1]].deref(i);
            }
        }
        if (sum < 0.99)printf("got rowsum of %e\n", sum);

    }
    tfor.put(daout, (unsigned char)0, (unsigned char)255);
    /*for(coor[1]=0;coor[1] < daouf.dims[1];coor[1]++){
        for(coor[0]=0;coor[0] < daouf.dims[0];coor[0]++){
            daouf(coor) = par_rowscale[coor[1]][coor[0]];
        }
    }
    tfos.put(daouf, (float)0, (float)255);*/


    coor[0] = par_RM.sizes[1];
    coor[1] = colhid.getNBcols();
    daout.setSizes(coor);
    coor[0] = 2;
    daouf.setSizes(coor);
    optr = daout.data;
    for(coor[1]=0;coor[1] < daout.dims[1];coor[1]++){
        if (colhid.data[coor[1]].getSize() == 0) printf("error! col %i has no state!\n", coor[1]);
        sum =0;
        for(coor[0]=0;coor[0] < daout.dims[0];coor[0]++){
            if ((i = colhid.data[coor[1]].find(coor[0])) == 0xFFFFFFFF) *(optr++) = 0;
            else {
                sum += colhid.data[coor[1]].deref(i);
                *(optr++) = (unsigned char)(colhid.data[coor[1]].deref(i) * 255);
            }
        }
         if (sum < 0.99) printf("got rowsum of %e\n", sum);
    }
    tfoc.put(daout, (unsigned char)0, (unsigned char)255);
    /*for(coor[1]=0;coor[1] < daouf.dims[1];coor[1]++){
        for(coor[0]=0;coor[0] < daouf.dims[0];coor[0]++){
            daouf(coor) = par_colscale[coor[1]][coor[0]];
        }
    }
    tfos.put(daouf, (float)0, (float)255);*/


    coor[0] = rowhid.getNBcols();
    coor[1] = colhid.getNBcols();
    daouf.setSizes(coor);

    double tumbufp[2];
    for(coor[1]=0;coor[1] < daouf.dims[1];coor[1]++){
        Tuple<Tuple<double,2u> > curproj = par_RM * colhid.getColumn(coor[1]);

        for(coor[0]=0;coor[0] < daouf.dims[0];coor[0]++){
            Tuple<double,2u> coeff;
            coeff = curproj.mkInnerProd(rowhid.getColumn(coor[0]) );
            coeff += par_rowscale[coor[0]];
            coeff += par_colscale[coor[1]];
            if (ExOp::isValid(coeff[0])){
                deflatedNB_getMeanAndVar(tumbufp, coeff[0],coeff[1]);
                daouf(coor) = tumbufp[0];
            }else daouf(coor) = 1.0;
        }
    }
    tfom.put(daouf, (float)0, (float)255);
    for(coor[1]=0;coor[1] < daouf.dims[1];coor[1]++){
        Tuple<Tuple<double,2u> > curproj = par_RM * colhid.getColumn(coor[1]);
        for(coor[0]=0;coor[0] < daouf.dims[0];coor[0]++){
            Tuple<double,2u> coeff;
            coeff = curproj.mkInnerProd(rowhid.getColumn(coor[0]) );
            coeff += par_rowscale[coor[0]];
            coeff += par_colscale[coor[1]];
            if (ExOp::isValid(coeff[0])) daouf(coor) = coeff[0];
            else daouf(coor) = 0.0;
        }
    }
    tfom.put(daouf, (float)0, (float)255);
    for(coor[1]=0;coor[1] < daouf.dims[1];coor[1]++){
        Tuple<Tuple<double,2u> > curproj = par_RM * colhid.getColumn(coor[1]);
        for(coor[0]=0;coor[0] < daouf.dims[0];coor[0]++){
            Tuple<double,2u> coeff;
            coeff = curproj.mkInnerProd(rowhid.getColumn(coor[0]) );
            coeff += par_rowscale[coor[0]];
            coeff += par_colscale[coor[1]];
            if (ExOp::isValid(coeff[0])) daouf(coor) = coeff[1];
            else daouf(coor) = 0.0;
        }
    }
    tfom.put(daouf, (float)0, (float)255);


}
LFHTEMP void BiClassifier<C,R,S>::exportZscore(const char* folderpath, const SparseMatrix<uint32_t> &input)const{
    DataGrid<float,2u> imout[7];
    Tuple<uint32_t,2u> coor;
    char buffer[1024];
    uint32_t i = strlen(folderpath);
    memcpy(buffer, folderpath, i);
    strcpy(buffer+i, "Zscore.tif");
    TiffFile tfoz(buffer,true);
    //coor3[0] = 3;
    Vector<uint32_t> rows = input.compute_rowsizes();
    coor[0] = rows.getSize();

    coor[1] = input.getNBcols();
    imout[0].setSizes(coor);
    imout[1].setSizes(coor);
    imout[2].setSizes(coor);
    imout[3].setSizes(coor);
    imout[4].setSizes(coor);
    imout[5].setSizes(coor);
    imout[6].setSizes(coor);
    uint32_t ite,k;
    double bufout[6];
    double tmp;
    double buf[2];
    for(coor[1] =0; coor[1]< imout[0].dims[1];coor[1]++){
        Tuple<Tuple<double,2u> > curproj = par_RM * colhid.getColumn(coor[1]);
        for(coor[0] =0; coor[0]<imout[0].dims[0];coor[0]++){
            Tuple<double,2u> coeff = curproj.mkInnerProd(rowhid.getColumn(coor[0]) );
            coeff += par_rowscale[coor[0]];
            coeff += par_colscale[coor[1]];
            if ((ite = input.data[coor[1]].find(coor[0])) == 0xFFFFFFFF) {imout[0](coor) =0.0f;
                if (ExOp::isValid(coeff[0])) deflatedNB_getMeanAndVar(buf, coeff[0],coeff[1]);
                else {buf[0] = 1.0; buf[1] =0.0;}
                imout[1](coor) =buf[0];imout[2](coor) =0.0f;imout[3](coor) =0.0f;imout[4](coor) =0.0f;imout[5](coor) =0.0f;imout[6](coor) =0.5f;
            }else{
                imout[0](coor) = input.data[coor[1]].deref(ite);
                if (ExOp::isValid(coeff[0])) deflatedNB_getMeanAndVar(buf, coeff[0],coeff[1]);
                else {buf[0] = 1.0; buf[1] =0.0;}
                // mean =  (1 - p(k = 0))

                imout[2](coor) = buf[0] - input.data[coor[1]].deref(ite);
                imout[3](coor) = log(buf[0]) - log((double)input.data[coor[1]].deref(ite));
                imout[1](coor) = buf[0];
                imout[4](coor) = (buf[0] - input.data[coor[1]].deref(ite)) / sqrt(buf[1]); //colhid.getColumn(coor[1])[0] * rowhid.getColumn(coor[0])[1] - log((double)input.data[coor[1]].deref(ite));
                if (ExOp::isValid(coeff[0])){
                    deflatedNBterm_wr_f_d2(bufout, input.data[coor[1]].deref(ite), coeff[0], coeff[1]);

                    imout[5](coor) = bufout[0] - lngamma(1.0 + input.data[coor[1]].deref(ite));
                    imout[6](coor) = deflatedNB_pvalue(input.data[coor[1]].deref(ite), coeff[0], coeff[1]);
                }else{
                    imout[5](coor) = 0.0;
                    imout[6](coor) = 0.0;
                }
            }
        }
    }

    tfoz.put(imout[0],(float)0.0f,(float)1.0f); // input
    tfoz.put(imout[1],(float)0.0f,(float)1.0f); // ???
    tfoz.put(imout[2],(float)0.0f,(float)1.0f); // difference to expectation
    tfoz.put(imout[3],(float)0.0f,(float)1.0f); // log difference to expectation
    tfoz.put(imout[4],(float)0.0f,(float)1.0f); // Z_score
    tfoz.put(imout[5],(float)0.0f,(float)1.0f); // Log likelihood of observation
    tfoz.put(imout[6],(float)0.0f,(float)1.0f); // p-value
}

LFHTEMP SparseMatrix<double> BiClassifier<C,R,S>::computeDeviation(const SparseMatrix<uint32_t> &input)const{ SparseMatrix<double> fout;
    fout.setNBcols(input.data.getSize());
    Tuple<uint32_t,2u> coor;
    uint32_t i;
    double buf[2];
    for(coor[1] = 0; coor[1]< input.data.getSize();coor[1]++){
        Tuple<Tuple<double,2u> > curproj = par_RM * colhid.getColumn(coor[1]);
        for(i =0; i<input.data[coor[1]].getSize();i++){
            coor[0] = input.data[coor[1]].deref_key(i);
            Tuple<double,2u> coeff;
            coeff = curproj.mkInnerProd(rowhid.getColumn(coor[0]) );
            coeff += par_colscale[coor[1]];
            coeff += par_rowscale[coor[0]];
            deflatedNB_getMeanAndVar(buf, coeff[0],coeff[1]);
            fout.data[coor[1]][coor[0]] = (((double)input.data[coor[1]].deref(i)) - buf[0]);
        }
    }
return fout;}

LFHTEMP class BiClassifier<C,R,S>::normalizeMK2_task{
    public:
        SparseMatrix<double> &output;
        const SparseMatrix<uint32_t> &input;
        const BiClassifier<C,R,S> &target;
        Tuple<uint32_t> colrange;
     normalizeMK2_task(const BiClassifier<C,R,S> &_target, const SparseMatrix<uint32_t> &_input, SparseMatrix<double> &_output, uint32_t nbthreads):target(_target),input(_input),output(_output){
        colrange.setSize(nbthreads+1);
        int i;
        for(i=0;i<nbthreads;i++) {
            colrange[i] = (i * input.getNBcols()) / nbthreads;
        }
        colrange[i] = input.getNBcols();
    }

    uint32_t operator()(uint32_t threadID){
        Tuple<Tuple<double,2u> > curproj;
        Tuple<double,2u> coeff;
        int i;
        uint32_t row;
        for(uint32_t col = colrange[threadID]; col < colrange[threadID+1]; col++){thrbase.updateProgress(threadID);
            if (auto ite =  target.colhid.getColumn(col).getIterator()) do {if (ite() >= target.par_RM.sizes[1]) printf("coc c%i!\n", ite());}while(ite++);
            curproj = target.par_RM * target.colhid.getColumn(col);
            for(i =0; i<input.data[col].getSize();i++){
                row = input.data[col].deref_key(i);
                if (row > target.rowhid.getNBcols()) printf("datarow is too far!\n");
                if (auto ite =  target.rowhid.getColumn(row).getIterator()) do {if (ite() >= target.par_RM.sizes[0]) printf("ror r%i!\n", ite());}while(ite++);
                coeff = curproj.mkInnerProd(target.rowhid.getColumn(row) );
                if (ExOp::isValid(target.par_rowscale[row])) coeff += target.par_rowscale[row];
                if (ExOp::isValid(target.par_colscale[col])) coeff += target.par_colscale[col];
                output.data[col][row] = LogPvalue_to_stdnorm(LogPvalue_NBdistrib_exppara_LH_exact(input.data[col].deref(i), coeff[0], -coeff[1]));
            }
        }
        thrbase.finishProgress(threadID);
    return 0;}

};
LFHTEMP SparseMatrix<double> BiClassifier<C,R,S>::normalizeMK2(const SparseMatrix<uint32_t> &input)const{SparseMatrix<double> fout;
     BiClassifier<C,R,S>::normalizeMK2_task lotask(*this, input, fout, thrbase.nbthreads);
     fout.setNBcols(input.data.getSize());

    thrbase.startProgress("Normalizing", lotask.colrange.last());
    for(int i=thrbase.nbthreads-1;i>0;i--) thrbase.submit(lotask , i);
    thrbase.submit_ThenWait(lotask, 0);
return fout;}


LFHTEMP SparseMatrix<double> BiClassifier<C,R,S>::normalize(const SparseMatrix<uint32_t> &input, bool renormalize_pvalue_cdf, SparseMatrix<double>* optout_pval)const{ SparseMatrix<double> fout;
    fout.setNBcols(input.data.getSize());
    if (optout_pval) optout_pval->setNBcols(input.data.getSize());
    KeyElem<double, Tuple<uint32_t,2u> > val_input;
    uint32_t i,k;

    // X ~ ZDNB(R,M)  -> D ~ N(0,1)    yes but... dropout?
    // X ~ NB(R,M)  -> D ~ N(0,1) but the 0 deviations are not stored, but deviation and bio-dropout probability can be recovered

    // hence was is that suffitient information to get deviation and weight

    double buf[6];
    HeapTree<KeyElem<double, Tuple<uint32_t, 2u> > > vals;
    if (renormalize_pvalue_cdf){
        for(val_input.d[1] = 0; val_input.d[1]< input.data.getSize();val_input.d[1]++){
            Tuple<Tuple<double,2u> > curproj = par_RM * colhid.getColumn(val_input.d[1]);
            for(i =0; i<input.data[val_input.d[1]].getSize();i++){
                val_input.d[0] = input.data[val_input.d[1]].deref_key(i);
                Tuple<double,2u> coeff;
                coeff = curproj.mkInnerProd(rowhid.getColumn(val_input.d[0]) );
                if (ExOp::isValid(par_colscale[val_input.d[1]][0])) coeff += par_colscale[val_input.d[1]];
                if (ExOp::isValid(par_rowscale[val_input.d[0]][0])) coeff += par_rowscale[val_input.d[0]];
                val_input.k = deflatedNB_pvalue(input.data[val_input.d[1]].deref(i), coeff[0], coeff[1]);
                if (optout_pval) optout_pval->data[val_input.d[1]][val_input.d[0]] = val_input.k;
                vals.insert(val_input);
            }
        }
        buf[1] = 1.0 / vals.getSize();
        buf[0] = buf[1] * 0.5;
        while(!vals.isEmpty()){
            val_input = vals.pop();
            fout.data[val_input.d[1]][val_input.d[0]] = Pvalue_to_stdnorm(buf[0]);
            buf[0] += buf[1];
        }
    }else{
        for(val_input.d[1] = 0; val_input.d[1]< input.data.getSize();val_input.d[1]++){
            Tuple<Tuple<double,2u> > curproj = par_RM * colhid.getColumn(val_input.d[1]);
            for(i =0; i<input.data[val_input.d[1]].getSize();i++){
                val_input.d[0] = input.data[val_input.d[1]].deref_key(i);
                Tuple<double,2u> coeff;
                coeff = curproj.mkInnerProd(rowhid.getColumn(val_input.d[0]) );
                if (ExOp::isValid(par_colscale[val_input.d[1]][0])) coeff += par_colscale[val_input.d[1]];
                if (ExOp::isValid(par_rowscale[val_input.d[0]][0])) coeff += par_rowscale[val_input.d[0]];
                val_input.k = deflatedNB_pvalue(input.data[val_input.d[1]].deref(i), coeff[0], coeff[1]);
                if (optout_pval) optout_pval->data[val_input.d[1]][val_input.d[0]] = val_input.k;
                fout.data[val_input.d[1]][val_input.d[0]] = Pvalue_to_stdnorm(val_input.k);
            }
        }
    }
return fout;}

LFHTEMP SparseMatrix<double> BiClassifier<C,R,S>::normalize_dropout(const SparseMatrix<uint32_t> &input, SparseMatrix<double>* optout_logitpvalue)const{ SparseMatrix<double> fout;
    fout.setNBcols(input.data.getSize());
    if (optout_logitpvalue) optout_logitpvalue->setNBcols(input.data.getSize());
    KeyElem<double, Tuple<uint32_t,2u> > val_input;
    uint32_t i,k;
    double buf[6];
    double dropfact;
    for(val_input.d[1] = 0; val_input.d[1]< input.data.getSize();val_input.d[1]++){
        Tuple<Tuple<double,2u> > curproj = par_RM * colhid.getColumn(val_input.d[1]);
        for(i =0; i<input.data[val_input.d[1]].getSize();i++){
            val_input.d[0] = input.data[val_input.d[1]].deref_key(i);
            Tuple<double,2u> coeff;
            coeff = curproj.mkInnerProd(rowhid.getColumn(val_input.d[0]) );
            if (ExOp::isValid(par_colscale[val_input.d[1]][0])) coeff += par_colscale[val_input.d[1]];
            if (ExOp::isValid(par_rowscale[val_input.d[0]][0])) coeff += par_rowscale[val_input.d[0]];
            val_input.k = NB_logitpvalue(input.data[val_input.d[1]].deref(i), coeff[0], coeff[1] , par_coldropterm[val_input.d[1]] + par_rowdropterm[val_input.d[0]] );
            if (optout_logitpvalue) optout_logitpvalue->data[val_input.d[1]][val_input.d[0]] = val_input.k;
            fout.data[val_input.d[1]][val_input.d[0]] = logitPval_to_Zscore(val_input.k);
        }
    }
return fout;}

/*
LFHTEMP void BiClassifier<C,R,S>::BiSparse_scope::initRowIndexFromHidden(const SparseMatrix<double, 0u> &rowhid){
    Tuple<uint32_t, 2u> coor;
    double m;
    uint32_t b;
    rowindex.setSize(rowhid.getNBcols());
    for(coor[1]=0;coor[1]< rowhid.dims[1];coor[1]++){
        coor[0]=0; b=0; m = rowhid(coor);
        for(coor[0]++;coor[0]< rowhid.dims[0];coor[0]++){
            if (rowhid(coor) > m) { b = coor[0]; m = rowhid(coor);}
        }
        rowindex[coor[1]] =b;
    }
}
LFHTEMP void BiClassifier<C,R,S>::BiSparse_scope::initColIndexFromHidden(const SparseMatrix<double, 0u> &colhid){
    Tuple<uint32_t, 2u> coor;
    double m;
    uint32_t b;
    colindex.setSize(colhid.getNBcols());
    for(coor[1]=0;coor[1]< colhid.dims[1];coor[1]++){
        coor[0]=0; b=0; m = colhid(coor);
        for(coor[0]++;coor[0]< colhid.dims[0];coor[0]++){
            if (colhid(coor) > m) { b = coor[0]; m = colhid(coor);}
        }
        colindex[coor[1]] =b;
    }
}*/


LFHTEMP BiClassifier<C,R,S>::Shared_EMSCOPE::Shared_EMSCOPE(BiClassifier<C,R,S>& _target, const SparseMatrix<C>& _data, const Vector<uint32_t> &excl_list): target(_target){
//LFHTEMP template<int F> void BiClassifier<C,R,S>::Shared_EMSCOPE::initCounts(const SparseMatrix<C, F>& data, SparseMatrix<double, 0u> &rowhid, SparseMatrix<double, 0u> &colhid){
    if (excl_list.getSize() == 0){
        data = _data;
    }else{
        data.setNBcols(_data.getNBcols());
        myHashmap<uint32_t, void> which;
        for(uint32_t i=0;i<excl_list.getSize();i++) which.addEntry(excl_list[i]);
        for(uint32_t i=0;i<_data.getNBcols();i++){
            if (auto ite = _data.data[i].getIterator()) do{
                if (which.find(ite()) == 0xFFFFFFFF) data.data[i][ite()] = *ite;
            } while(ite++);
        }
    }
    trdata = data.mkTranspose();
    printf("fixing %i and %i nbrows so they match\n", trdata.data.getSize(), target.rowhid.getNBcols());
    while(trdata.data.getSize() < _target.rowhid.getNBcols()) trdata.data.push_back();

    init_routine(_target);
}
LFHTEMP BiClassifier<C,R,S>::Shared_EMSCOPE::Shared_EMSCOPE(BiClassifier<C,R,S>& _target, const SparseMatrix<Tuple<uint32_t, 2u> >& _data, const Vector<uint32_t> &excl_list, bool do_incl_introns): target(_target){
//LFHTEMP template<int F> void BiClassifier<C,R,S>::Shared_EMSCOPE::initCounts(const SparseMatrix<C, F>& data, SparseMatrix<double, 0u> &rowhid, SparseMatrix<double, 0u> &colhid){
    data.setNBcols(_data.getNBcols());
    introns.setNBcols(_data.getNBcols());
    myHashmap<uint32_t, void> which;
    for(uint32_t i=0;i<excl_list.getSize();i++) which.addEntry(excl_list[i]);

    if (do_incl_introns){
        for(uint32_t i = 0 ; i< _data.getNBcols();i++){
            if (auto ite = _data.data[i].getIterator()) do{
                if (which.find(ite()) != 0xFFFFFFFF) continue;
                data.data[i][ite()] = (*ite)[0] + (*ite)[1];
                introns.data[i][ite()] = (*ite)[1];
            }while(ite++);
        }
    }else{
        for(uint32_t i = 0 ; i< _data.getNBcols();i++){
            if (auto ite = _data.data[i].getIterator()) do{
                if (which.find(ite()) != 0xFFFFFFFF) continue;
                data.data[i][ite()] = (*ite)[0];
                introns.data[i][ite()] = (*ite)[1];
            }while(ite++);
        }
    }
    trdata = data.mkTranspose();
    printf("fixing %i and %i nbrows so they match\n", trdata.data.getSize(), target.rowhid.getNBcols());
    while(trdata.data.getSize() < _target.rowhid.getNBcols()) trdata.data.push_back();
    init_routine(_target);
}
LFHTEMP void BiClassifier<C,R,S>::Shared_EMSCOPE::init_routine(BiClassifier<C,R,S>& _target){
    Tuple<uint32_t, 2u> coor;
    ExOp::show(_target.par_RM.sizes);
    printf("das %i %i\n", data.getNBcols(), trdata.getNBcols()); fflush(stdout);
    coor[0] = _target.par_RM.sizes[0];
    coor[1] = data.getNBcols();
    col_nbmax_hrow.toSizes(coor).toZero();

    coor[0] = _target.par_RM.sizes[1];
    coor[1] = trdata.getNBcols();
    row_nbmax_hcol.toSizes(coor).toZero();

    coor[0] = _target.par_RM.sizes[0];
    coor[1] = _target.par_RM.sizes[1];
    maxgrid_nonzero_count.setSizes(coor).toZero();

    nb_maxrowh.setSize(_target.par_RM.sizes[0]).toZero();
    nb_maxcolh.setSize(_target.par_RM.sizes[1]).toZero();


    center_shift.setSize( (_target.par_RM.sizes[0]) > _target.par_RM.sizes[1] ? _target.par_RM.sizes[0] : _target.par_RM.sizes[1] ).toZero();

    istimefor_columns = true;
//    printf("Successfully allocated! %i\n", _target.rowhid.getNBcols()); fflush(stdout);
    uint32_t i,j;
    for(coor[1]=0;coor[1]<_target.rowhid.getNBcols();coor[1]++) {
        if (_target.rowhid.data[coor[1]].getSize() == 0) {printf("nothing for row %i\n", coor[1]);continue;}
        _target.rowhid.makeMaximumAsFirst(coor[1]);
        nb_maxrowh[_target.rowhid.data[coor[1]].deref_key(0)]++;
    }
//    printf("Successfully A! %i\n", _target.colhid.getNBcols()); fflush(stdout);
    for(coor[1]=0;coor[1]<_target.colhid.getNBcols();coor[1]++){
        if (_target.colhid.data[coor[1]].getSize() == 0) {printf("nothing for col %i\n", coor[1]);continue;}
        _target.colhid.makeMaximumAsFirst(coor[1]);
        nb_maxcolh[_target.colhid.data[coor[1]].deref_key(0)]++;
    }
//    printf("Successfully B!\n"); fflush(stdout);

    if (auto ite = data.getIterator()) do{
        if (_target.colhid.data[ite.getCol()].getSize() == 0) printf("impossible... col?\n");
        if (_target.rowhid.data[ite.getRow()].getSize() == 0) printf("impossible... row?\n");
        coor[0] = _target.colhid.data[ite.getCol()].deref_key(0);
        coor[1] = ite.getRow();
        row_nbmax_hcol(coor)++;

        coor[0] = _target.rowhid.data[ite.getRow()].deref_key(0);
        coor[1] = ite.getCol();
        col_nbmax_hrow(coor)++;

        coor[1] = _target.colhid.data[ite.getCol()].deref_key(0);
        maxgrid_nonzero_count(coor)++;

    }while (ite++);
//    printf("Successfully Exit!\n"); fflush(stdout);
}

/*
LFHTEMP template<int F> void BiClassifier<C,R,S>::BiSparse_scope::changeRowIndex(const SparseMatrix<C, F>& data, uint32_t rowId, uint32_t value){
    Tuple<uint32_t, 2u> coor;
    for(coor[1]=0;coor[1]<data.getNBcols();coor[1]++){
        if (data.data[coor[1]].find(rowId) != 0xFFFFFFFF){
            coor[0] = rowindex[rowId]; rowh_per_coldata(coor)--;
            coor[0] = value; rowh_per_coldata(coor)++;
        }
    }
    rowindex[rowId] = value;
}
LFHTEMP template<int F> void BiClassifier<C,R,S>::BiSparse_scope::changeColIndex(const SparseMatrix<C, F>& data, uint32_t colId, uint32_t value){
    Tuple<uint32_t, 2u> coor;
    int32_t i;
    for(i=0;i<data[colId].getSize();i++){
        coor[1] = data[colId].deref_key(i);
        coor[0] = colindex[colId]; colh_per_rowdata(coor)--;
        coor[0] = value; colh_per_rowdata(coor)++;
    }
    colindex[colId] = value;
}*/
LFHTEMP void BiClassifier<C,R,S>::Parallel_EMSCOPE::setRange(uint32_t r_start, uint32_t r_max,uint32_t c_start, uint32_t c_max){
    colrange[0] = c_start; colrange[1]=c_max;
    rowrange[0] = r_start; rowrange[1]=r_max;
}
/*LFHTEMP void BiClassifier<C,R,S>::Parallel_EMSCOPE::setPointers(Tuple<double>* _r, Tuple<double>* _c , uint32_t* _h){
    rowdata=_r;coldata=_c;cachind =_h;
}*/
LFHTEMP void BiClassifier<C,R,S>::Parallel_EMSCOPE::setSizes(Tuple<uint32_t, 2u> &coor, int nbrow){
    dll_dMR2.setSizes(coor[0],coor[1]);
    Tuple<uint32_t, 2u> ncoor;
    ncoor[0] = coor[0];
    ncoor[1] = coor[1];
    nonzero_diff.setSizes(coor);
    nbmax_diff.setSize(coor[0] > coor[1]? coor[0] : coor[1]);
    dd_massive.setSize(coor[1]*coor[0]*2);
}
LFHTEMP uint32_t BiClassifier<C,R,S>::Parallel_EMSCOPE::operator()(){
    ExOp::toZero(ll);
    dd_massive.toZero();
    dll_dMR2.toZero();
    nonzero_diff.toZero();
    nbmax_diff.toZero();
    indexes_changes.toMemfree();
    Tuple<uint32_t,2u> coor, acoor;

    uint32_t i,j;
    Tuple<Tuple<double,2u> > curproj[4]; // R H_j C_j
    double tmp;
    Tuple< Tuple<double, 2u>, 0u> future_statistics[3];
    Trianglix< Tuple<double, 3u>, 0u> future_statistic_trix[3];
    Tuple<double, 3u> trix_input;
    SparseTuple< Tuple<double, 2u> > future_factor;

    SparseTuple<double> candidates[2];
    SparseTuple<double>  can_backup;
    Tuple<double, 3u> lls;
    Tuple<double, 3u> lld;
    Tuple<double, 2u> coeff;
    Tuple<double, 2u> coeff_test;
    double sum[4];

    FunctionScope<double&,double&,double> funabsoft(ExCo<double>::toAbsoftInverse,pow(0.5, 256.0f));

    Tuple<double, 2u> d_scale, o_scale, d_scale_tmp;
    Trianglix<double, 2u> dd_scale;

    Trianglix<double, 2u> dd_scale_tmp;

    Tuple<Tuple<double,2u>, 0u> rever;

    Tuple<double> prior_deriv;
    double prior_llterm[3];
    double Rval,Mval;

    double logexpmin;
    double expRterm;
    union LocalDef{
        double evilterms[6];
        KeyElem<double, Tuple<double, 5u> > evilterm_alt;
        LocalDef(){}
    };
    LocalDef ld;

    scp->target.par_RM.getDims(coor);
    Tuple<Tuple<double,2u> > RM_deriv;

    Tuple<Tuple<Tuple<double,2u> > > rever_buf;
    Tuple<Tuple<double, 6u> > dd_buf;

    uint32_t dir;
    bool haschange;

    int insist=0;

    double prediction;

    double accuracy =0.0;
    nb_fail =0u;
    nb_insist =0u;
    uint32_t nb_funky=0u;

    Tuple<double, 4u> tosolve_tmp[3];

    Tuple<double, 3u> tosolve_num2;
    Tuple<double, 3u> tosolve_num2_tmp;
    Trianglix<double, 3u> tosolve_den2;
    Trianglix<double, 3u> tosolve_den2_tmp;
    Tuple<double> drop_buf;

    double expexp;
    double expexp2;
    bool checksafe;
    Vector<int> torem;
    uint32_t toremdirection;
    ProgressBarPrint pbar(20);
    if (use_stdout) pbar.start((scp->istimefor_columns) ? "Compute columns": "Compute rows");

    if (scp->istimefor_columns){

        RM_deriv.setSize(scp->target.par_RM.sizes[1]);
        prior_deriv.setSize(scp->target.par_RM.sizes[1]);

        rever_buf.setSize(scp->target.rowhid.getNBcols());
        dd_buf.setSize(scp->target.rowhid.getNBcols());

        for(i=0;i<3;i++) {
            future_statistics[i].setSize(scp->target.par_RM.sizes[0]);
            future_statistic_trix[i].setSize(scp->target.par_RM.sizes[0]);
        }
        ExOp::toZero(sum);
        for(coor[1]=colrange[0];coor[1] < colrange[1];coor[1]++){
            if (use_stdout) pbar.update(((double)coor[1] - colrange[0])/ (colrange[1]-colrange[0]));
            //printf("cur coor %i/%i\n", coor[1], scp->target.colhid.getNBcols());fflush(stdout);
            const SparseTuple<double>& curcol = scp->target.colhid.getColumn(coor[1]);
            if (curcol.getSize() == 0) {printf("Queried column %i is empty... impossible!\n",coor[1]); exit(1);}
            if (!ExOp::isValid(scp->target.par_colscale[coor[1]][0])){
                scp->newscale[coor[1]] = scp->target.par_colscale[coor[1]];

                if (scp->data.getColSize(coor[1]) == 0) {candidates[0] = scp->target.colhid.getColumn(coor[1]); scp->newhidden.memmoveColumn(coor[1], candidates[0]); continue;} // nothing to see...
                // column only has 1s, cannot be used to update R or M, drop out only, pure state guarrantied

                // put in stone... for now
                candidates[0] = scp->target.colhid.getColumn(coor[1]); scp->newhidden.memmoveColumn(coor[1], candidates[0]); continue;
                // put in stone... for now

                RM_deriv.toZero();

                for(coor[0]=0;coor[0]<scp->col_nbmax_hrow.dims[0];coor[0]++){
                    acoor[0] = coor[0];
                    for(acoor[1] =0;acoor[1]< scp->target.par_RM.sizes[1]; acoor[1]++){
                        RM_deriv[acoor[1]][0] += scp->col_nbmax_hrow(coor) * log(scp->target.par_D(acoor));
                        RM_deriv[acoor[1]][0] += (scp->nb_maxrowh[coor[0]] - scp->col_nbmax_hrow(coor)) * log(1.0 - scp->target.par_D(acoor));
                    }
                }
                acoor[1] =0;
                for(acoor[0]=1 ;acoor[0]< scp->target.par_RM.sizes[1]; acoor[0]++){
                    if (RM_deriv[acoor[0]][0] > RM_deriv[acoor[1]][0]) acoor[1] = acoor[0];
                }
                candidates[0].toZero();
                candidates[0][acoor[1]] = 1.0;
                if (candidates[0].deref_key(0) != curcol.deref_key(0)){
                    nbmax_diff[candidates[0].deref_key(0)]++;
                    nbmax_diff[curcol.deref_key(0)]--;
                    for(i = 0; i < scp->data.data[coor[1]].getSize();i++){
                        acoor[0] = scp->target.rowhid.data[ scp->data.data[coor[1]].deref_key(i)].deref_key(0);
                        acoor[1] = curcol.deref_key(0);
                        nonzero_diff(acoor)--;
                        acoor[1] = candidates[0].deref_key(0);
                        nonzero_diff(acoor)++;
                    }
                }
                scp->newhidden.memmoveColumn(coor[1], candidates[0]);
            }else{
                checksafe = false;

                prior_deriv = scp->target.colhid_prior * curcol;

                if (coor[1] >= scp->target.colhid.data.getSize()) exit(1);
                // Step 1 cache line specific data
                curproj[0] = scp->target.par_RM * curcol;

                //dll_dH.toZero();
                //dll_dC.toZero();
                // Step 2 find direction, find 2 vector for R and M
                RM_deriv.toZero();



                prior_llterm[0] = scp->target.colhid_prior.mkInnerMult(curcol) + scp->target.colhid_prior_constant;
                lls[2] = log(prior_llterm[0]) * ((scp->tmptmp) ? scp->data.data[coor[1]].getSize() : 1);
                //lls[2] = 0.0;
                future_statistics[2].toZero(); future_statistic_trix[2].toZero();
                d_scale.toZero(); dd_scale.toZero();
                for(i=0;i<scp->data.data[coor[1]].getSize();i++){
                    coor[0]= scp->data.data[coor[1]].deref_key(i);
                    if (!ExOp::isValid(scp->target.par_rowscale[coor[0]])) continue;
                    coeff = curproj[0].mkInnerProd(scp->target.rowhid.getColumn(coor[0]));
                    coeff += scp->target.par_rowscale[coor[0]];
                    coeff += scp->target.par_colscale[coor[1]];
                    deflatedNBterm_wr_f_d2(ld.evilterms, scp->data.data[coor[1]].deref(i)  , coeff[0], coeff[1]);
                    dd_scale.data[0] += ld.evilterms[3]; dd_scale.data[1] += ld.evilterms[5]; dd_scale.data[2] += ld.evilterms[4];
                    d_scale[0] += ld.evilterms[1]; d_scale[1] += ld.evilterms[2];

                    lls[2] += ld.evilterms[0];
                    //printf("(%e,%e), %e,%e,%e,%e,%e\n", scp->target.par_rowscale[coor[0]][0], scp->target.par_rowscale[coor[0]][1], ld.evilterms[0],ld.evilterms[1],ld.evilterms[2],ld.evilterms[3],ld.evilterms[4] );

                    //d_scale[0] += e += ld.evilterm_alt;
                    rever = scp->target.par_RM.mkBackMult(scp->target.rowhid.getColumn(coor[0]));

                    if (auto ite = rever.getIterator()) do{
                        RM_deriv[ite()][0] += (*ite)[0] * ld.evilterms[1];
                        RM_deriv[ite()][1] += (*ite)[1] * ld.evilterms[2];
                    }while(ite++);

                    if (auto ite = scp->target.rowhid.getColumn(coor[0]).getIterator()) do{

                        future_statistics[2][ite()][0] += ld.evilterm_alt.d[0] * (*ite);
                        future_statistics[2][ite()][1] += ld.evilterm_alt.d[1] * (*ite);
                        trix_input[0] = ld.evilterms[3] * (*ite);
                        trix_input[1] = ld.evilterms[4] * (*ite);
                        trix_input[2] = ld.evilterms[5] * (*ite) * 0.5;

                        if (auto ite2 = scp->target.rowhid.getColumn(coor[0]).getIterator()) do{
                            if (ite2() > ite()) continue;
                            future_statistic_trix[2].data[ite2() + ((ite() * (ite() + 1)) >> 1)] += trix_input * (*ite2);
                        }while(ite2++);
                    }while(ite++);

                    rever_buf[i].toMemmove(rever);
                    dd_buf[i][0] = ld.evilterms[0];
                    dd_buf[i][1] = ld.evilterms[1];
                    dd_buf[i][2] = ld.evilterms[2];
                    dd_buf[i][3] = ld.evilterms[3];
                    dd_buf[i][4] = ld.evilterms[4];
                    dd_buf[i][5] = ld.evilterms[5];

                } // dG/dz = c1f1 VM  dH/dz = c2 f2VR
                //printf("cur coor %i/%i\n", coor[1], scp->target.colhid.getNBcols());fflush(stdout);

                ll[0] += lls[2];
                prediction = 0.0f;
                // Step 3 (line search) define set of candidates


                sum[0] = 0.0;
                prior_llterm[2] = 0.0;
                candidates[0].toMemfree();
                for(i=0;i< RM_deriv.getSize();i++){
                    prior_llterm[1] = RM_deriv[i][0] + RM_deriv[i][1] + prior_deriv[i] * (2.0 * ((scp->tmptmp) ? scp->data.data[coor[1]].getSize() : 1) / prior_llterm[0]);
                    if ((curcol.find(i) != 0xFFFFFFFF)||(prior_llterm[1] > 0.0f)) {
                        candidates[0][i] = prior_llterm[1];
                        sum[0] +=  prior_llterm[1];
                    }
                }
                if (candidates[0].getSize() < 2u) { // stuck in a corner,... that's good but...
                    candidates[0] = curcol;
                    dd_scale_tmp = dd_scale; d_scale_tmp = d_scale;
                    d_scale = dd_scale.toEigenTransform(funabsoft) * d_scale;
                    d_scale *= scp->hiddenalpha;
                }else{
                    // project onto plane
                    coor[0]=0;
                    can_backup = candidates[0];
                    do{
                        torem.toMemfree();
                        sum[0] /= candidates[0].getSize();
                        if (auto ite = candidates[0].getIterator()) do{
                            if (((*ite) > sum[0])||(curcol.find(ite()) != 0xFFFFFFFF)) (*ite) -= sum[0];
                            else torem.push_back(ite());
                        }while(ite++);
                        if (torem.getSize() == 0u) {
                            sum[0] = 0.0;
                            if (auto ite = candidates[0].getIterator()) do{sum[0] += (*ite);}while(ite++);
                            if (fabs(sum[0]) < 0.0001) break;
                            if ((coor[0]++) == 5) {
                                printf("warning: got stuck... catastrophy cancellations\n");
                                can_backup.show();
                                curcol.show();
                                coor[0] = 0;
                                if (auto ite = candidates[0].getIterator()) do{if ((*ite)>= 0) coor[0]++;}while(ite++);
                                if (auto ite = candidates[0].getIterator()) do{(*ite) = ((*ite)>= 0) ?  candidates[0].getSize() - coor[0]: coor[0];}while(ite++);
                                if (auto ite = candidates[0].getIterator()) do{(*ite) /=  candidates[0].getSize();}while(ite++);
                                printf("warning: got stuck... catastrophy cancellations\n");
                                exit(1);
                                break;
                            }
                        }else{
                            if (candidates[0].getSize() - torem.getSize() <= 1){ // would leave nothing... try the highest second?
                                if (candidates[0].getSize() == torem.getSize()) {printf("Impossibleble!\n"); exit(1);}
                                acoor[0] = 0;
                                for(acoor[1]=1;acoor[1] < torem.getSize();acoor[1]++) if (candidates[0][torem[acoor[0]]] < candidates[0][torem[acoor[1]]]) acoor[0] = acoor[1];
                                for(i=0;i<torem.getSize();i++) {
                                    if (i == acoor[0]) continue;
                                    j = candidates[0].find(torem[i]);
                                    if (j == 0xFFFFFFFF) {printf("erasing unexisting... !\n"); exit(1);}
                                    candidates[0].erase_from_iterator(j);
                                }
                                 if (candidates[0].getSize() != 2) {printf("2 should remain. got %i\n", candidates[0].getSize()); exit(1);}
                                if (candidates[0].deref_key(0) == torem[acoor[0]]){
                                    candidates[0].deref(0) = 1.0;
                                    candidates[0].deref(1) = -1.0;
                                }else{
                                    candidates[0].deref(0) = -1.0;
                                    candidates[0].deref(1) = 1.0;
                                }
                                break;
                            }
                            for(i=0;i<torem.getSize();i++) {
                                j = candidates[0].find(torem[i]);
                                if (j == 0xFFFFFFFF) {printf("trying to remove un-existing!\n"); exit(1);}
                                candidates[0].erase_from_iterator(j);
                            }
                            sum[0] =0.0;
                            if (auto ite = candidates[0].getIterator()) do{sum[0] += *ite;}while(ite++);
                        }
                    }while(true);

                    if (auto ite = candidates[0].getIterator()){
                        while((*ite) > 0.0) ite++;
                        sum[0] = curcol[dir = ite()] / (*ite);

                        while(ite++){
                            if (*ite > 0.0) continue;
                            sum[1] = curcol[ite()] / (*ite);
                            if (sum[1] > sum[0]) {sum[0] = sum[1]; dir = ite();}
                        }
                        sum[0] = -sum[0]; // w has to be in 0-sum[0] range

                        tosolve_num2[0] = (scp->target.colhid_prior.mkInnerMult(candidates[0],curcol) + scp->target.colhid_prior_constant) * ( 2.0 * ((scp->tmptmp) ? scp->data.data[coor[1]].getSize() : 1)/ prior_llterm[0]);
                        tosolve_num2[1] = d_scale[0];
                        tosolve_num2[2] = d_scale[1];
                        tosolve_den2.data[1] = 0.0; tosolve_den2.data[3] = 0.0;
                        tosolve_den2.data[2] = dd_scale.data[0];
                        tosolve_den2.data[5] = dd_scale.data[2];
                        tosolve_den2.data[4] = dd_scale.data[1];
                        tosolve_den2.data[0] = (scp->target.colhid_prior.mkInnerMult(candidates[0])+ scp->target.colhid_prior_constant) * ( 2.0 * ((scp->tmptmp) ? scp->data.data[coor[1]].getSize() : 1) / prior_llterm[0]);
                        for(i=0;i<scp->data.data[coor[1]].getSize();i++){
                            coor[0]= scp->data.data[coor[1]].deref_key(i);
                            if (!ExOp::isValid(scp->target.par_rowscale[coor[0]])) continue;
                            coeff = rever_buf[i].mkInnerProd(candidates[0]);
                            tosolve_num2[0] += coeff[0] * dd_buf[i][1]+ coeff[1] * dd_buf[i][2];
                            tosolve_den2.data[0] += dd_buf[i][3] * coeff[0] * coeff[0] + dd_buf[i][4] * coeff[1] * coeff[1] + dd_buf[i][5] * coeff[0] * coeff[1];
                            tosolve_den2.data[1] += dd_buf[i][3] * coeff[0] + dd_buf[i][5] * coeff[1] * 0.5;
                            tosolve_den2.data[3] += dd_buf[i][4] * coeff[1] + dd_buf[i][5] * coeff[0] * 0.5;
                        }

                        tosolve_den2_tmp = tosolve_den2;
                        tosolve_num2_tmp = tosolve_num2;
                        tosolve_den2.toEigenTransform(funabsoft);
                        tosolve_num2 = tosolve_den2 * tosolve_num2;

                        tosolve_num2 *= scp->hiddenalpha;

                        expexp = tosolve_num2[0] * tosolve_num2_tmp[0] + tosolve_num2[1] * tosolve_num2_tmp[1] + tosolve_num2[2] * tosolve_num2_tmp[2];
                        expexp2 = tosolve_den2_tmp.data[0] * tosolve_num2[0] * tosolve_num2[0] * 0.5;
                        expexp2 += tosolve_den2_tmp.data[2] * tosolve_num2[1] * tosolve_num2[1] * 0.5;
                        expexp2 += tosolve_den2_tmp.data[5] * tosolve_num2[2] * tosolve_num2[2] * 0.5;
                        expexp2 += tosolve_den2_tmp.data[1] * tosolve_num2[0] * tosolve_num2[1];
                        expexp2 += tosolve_den2_tmp.data[3] * tosolve_num2[0] * tosolve_num2[2];
                        expexp2 += tosolve_den2_tmp.data[4] * tosolve_num2[1] * tosolve_num2[2];

                        if (tosolve_num2[0] <= 0.0){ // went the other way... just update scale
                            if (tosolve_num2[0] == 0.0){
                                candidates[0] = curcol;
                            }else{
                                candidates[0] = curcol;
                                dd_scale_tmp = dd_scale; d_scale_tmp = d_scale;
                                d_scale = dd_scale.toEigenTransform(funabsoft) * d_scale;
                                d_scale *= scp->hiddenalpha;

                            }
                       }else if (tosolve_num2[0] >= sum[0]){ // too far!, find solution on boundary...
                            d_scale[0] = tosolve_num2_tmp[1] + sum[0] * tosolve_den2.data[1];
                            d_scale[1] = tosolve_num2_tmp[2] + sum[0] * tosolve_den2.data[3];
                            d_scale = dd_scale.toEigenTransform(funabsoft) * d_scale;
                            d_scale *= scp->hiddenalpha;

                            candidates[1] = candidates[0];
                            if (auto ite = candidates[0].getIterator()) do{
                                (*ite) = (*ite) * sum[0]  + (((j = curcol.find(ite())) == 0xFFFFFFFF) ? 0.0: curcol.deref(j));
                            }while(ite++);

                            if (auto ite = curcol.getIterator()) do{
                                if (candidates[0].find(ite())== 0xFFFFFFFF) candidates[0][ite()] = *ite;
                            }while(ite++);
                            j = candidates[0].find(dir);
                            if (j == 0xFFFFFFFF) {printf("trying to remove un-existing!\n"); exit(1);}
                            candidates[0].erase_from_iterator(j);
                            sum[1] =0.0;
                            if (auto ite = candidates[0].getIterator()) {
                                do{
                                    sum[1] += *ite;
                                }while(ite++);
                                if (fabs(sum[1] - 1.0) > 0.1) {printf("bound proj dir %i... hidden state norm is %e! %e is fact\n", dir, sum[1], sum[0]);
                                candidates[1].show();
                                candidates[0].show();
                                curcol.show();
                                exit(1);}
                            }else {printf("got empty candidate!\n"); exit(1);}
                            candidates[1].toMemfree();
                        }else{ // all set!
                            //checksafe = true;

                            d_scale[0] = tosolve_num2[1];
                            d_scale[1] = tosolve_num2[2];
                            candidates[1] = candidates[0];
                            if (auto ite = candidates[0].getIterator()) do{
                                (*ite) = (*ite) * tosolve_num2[0] + (((j = curcol.find(ite())) == 0xFFFFFFFF) ? 0.0: curcol.deref(j));
                            }while(ite++);
                            if (auto ite = curcol.getIterator()) do{
                                if (candidates[0].find(ite())== 0xFFFFFFFF) candidates[0][ite()] = *ite;
                            }while(ite++);
                            sum[1] = 0.0;
                            if (auto ite = candidates[0].getIterator()) {
                                do{
                                    sum[1] += *ite;
                                }while(ite++);
                                if (fabs(sum[1] - 1.0) > 0.1) {printf("normal proj... hidden state norm is %e! %e is fact\n", sum[1], tosolve_num2[0]);
                                candidates[0].show();
                                candidates[1].show();
                                curcol.show();
                                LFH_ALIVE;exit(1);}
                            }else {printf("got empty candidate!\n"); exit(1);}
                            candidates[1].toMemfree();
                        }
                        candidates[0].makeMaximumAsFirst();
                        if (!scp->ignore_drop){
                            if (candidates[0].deref_key(0) != curcol.deref_key(0)) {
                                // major change, add *last* boundary on the search space!
                                tosolve_num2[2] = ((i = curcol.find(candidates[0].deref_key(0))) == 0xFFFFFFFF) ? 0.0 : curcol.deref(i);
                                if (auto ite = curcol.getIterator()) {
                                    tosolve_num2[0] = (*ite) - tosolve_num2[2]; // guarranied to be positive
                                    tosolve_num2[1] = candidates[0].deref(0) -  (((i = candidates[0].find(ite())) == 0xFFFFFFFF) ? 0.0 : candidates[0].deref(i));
                                    sum[1] = tosolve_num2[1] / (tosolve_num2[0] + tosolve_num2[1]);
                                    j = 0;
                                    while(ite++){
                                        if ((tosolve_num2[0] =  (*ite) - tosolve_num2[2]) <= 0) continue;
                                        tosolve_num2[1] = candidates[0].deref(0) -  (((i = candidates[0].find(ite())) == 0xFFFFFFFF) ? 0.0 : candidates[0].deref(i));
                                        sum[0] = tosolve_num2[1] / (tosolve_num2[0] + tosolve_num2[1]);
                                        if (sum[0] < sum[1]) {j = ite.getOffset(); sum[1] = sum[0];}
                                    }
                                }
                                sum[0] = candidates[0].deref(0) * (1.0 - sum[1]) + sum[1] * (((i = curcol.find(candidates[0].deref_key(0))) == 0xFFFFFFFF) ? 0.0 : curcol.deref(i));
                                candidates[1][curcol.deref_key(j)] = sum[0];
                                candidates[1][candidates[0].deref_key(0)] = sum[0];
                                if (auto ite = curcol.getIterator()) do {
                                    if (candidates[1].find(ite()) == 0xFFFFFFFF) candidates[1][ite()] = (*ite) * sum[1] + (1.0 - sum[1]) * (((i = candidates[0].find(ite())) == 0xFFFFFFFF) ? 0.0 : candidates[0].deref(i));
                                }while(ite++);
                                if (auto ite = candidates[0].getIterator()) do {
                                    if (candidates[1].find(ite()) == 0xFFFFFFFF) candidates[1][ite()] = (*ite) * (1.0 - sum[1]);
                                }while(ite++);

                                if (auto ite = candidates[1].getIterator()) {
                                    sum[0] = 0.0;
                                    do{
                                        sum[0] += *ite;
                                    }while(ite++);
                                    if (fabs(sum[0] - 1.0) > 0.1) {printf("hidden state norm cand1 is %e!\n", sum[0]);
                                        candidates[1].show();
                                        printf("frac %e\n", sum[1]);
                                        candidates[0].show();
                                        curcol.show();
                                    LFH_ALIVE;exit(1);}
                                }
                            }/*else{ // forcefully try closest maybe?

                            }*/
                        }
                    }else{
                        // stuck in corner again...
                        candidates[0] = curcol;
                        dd_scale_tmp = dd_scale; d_scale_tmp = d_scale;
                        d_scale = dd_scale.toEigenTransform(funabsoft) * d_scale;
                        d_scale *= scp->hiddenalpha;

                    }
                }
                sum[0] = 0.0f;
                if (auto ite = candidates[0].getIterator()) {
                    do{
                        sum[0] += *ite;
                    }while(ite++);
                    if (fabs(sum[0] - 1.0) > 0.1) {printf("hidded state norm is hh %e!\n", sum[0]);
                        candidates[0].show();
                        candidates[1].show();
                    LFH_ALIVE;exit(1);}
                }else {printf("got empty candidate!\n"); exit(1);}

                // find intersection on line
                // new scale!
                //printf("cur coor %i/%i\n", coor[1], scp->target.colhid.getNBcols());fflush(stdout);

                coeff = scp->target.par_colscale[coor[1]];
                coeff += d_scale;
                if ((tmp = (coeff[0]*coeff[0] + coeff[1]*coeff[1])) > 100.0){ // too far, send to boundary
                    tmp = pow(tmp, -0.5);
                    d_scale[0] = coeff[0] * tmp;
                    d_scale[1] = coeff[1] * tmp;
                    d_scale -= scp->target.par_colscale[coor[1]];
                }

                // Step 4 Evaluate Quantidates

                insist =0;
                while(true){

                    for(j=0;j<2;j++){
                        future_statistics[j].toZero(); future_statistic_trix[j].toZero();
                        curproj[j] = scp->target.par_RM * candidates[j];
                        lls[j] = log(scp->target.colhid_prior.mkInnerMult(candidates[j]) + scp->target.colhid_prior_constant) * ((scp->tmptmp) ? scp->data.data[coor[1]].getSize() : 1);
                       // lls[j] = 0.0;
                        if ((candidates[1].getSize() == 0)||(insist > 0)) break;
                    }
                    for(i=0;i<scp->data.data[coor[1]].getSize();i++){
                        coor[0]= scp->data.data[coor[1]].deref_key(i);
                        if (!ExOp::isValid(scp->target.par_rowscale[coor[0]])) continue;
                        for(j=0;j<2;j++){
                            coeff = curproj[j].mkInnerProd(scp->target.rowhid.getColumn(coor[0]));
                            coeff += scp->target.par_rowscale[coor[0]];
                            coeff += scp->target.par_colscale[coor[1]];
                            coeff += d_scale;


                            deflatedNBterm_wr_f_d2(ld.evilterms, scp->data.data[coor[1]].deref(i), coeff[0], coeff[1]);
                            if (auto ite = scp->target.rowhid.getColumn(coor[0]).getIterator()) do{
                                future_statistics[j][ite()][0] += ld.evilterm_alt.d[0] * (*ite);
                                future_statistics[j][ite()][1] += ld.evilterm_alt.d[1] * (*ite);
                                trix_input[0] = ld.evilterms[3] * (*ite);
                                trix_input[1] = ld.evilterms[4] * (*ite);
                                trix_input[2] = ld.evilterms[5] * (*ite) * 0.5;
                                if (auto ite2 = scp->target.rowhid.getColumn(coor[0]).getIterator()) do{
                                    if (ite2() > ite()) continue;
                                    future_statistic_trix[j].data[ite2() + ((ite() * (ite() + 1)) >> 1)] += trix_input * (*ite2);
                                }while(ite2++);
                            }while(ite++);

                            lls[j] += ld.evilterms[0];
                            if ((candidates[1].getSize() == 0)||(insist > 0)) break;
                        }
                    } // H(z) = fcfc vMMv

                    if (checksafe){
                        double diffdiff = (log(scp->target.colhid_prior.mkInnerMult(candidates[0])+ scp->target.colhid_prior_constant) - log(scp->target.colhid_prior.mkInnerMult(curcol)+ scp->target.colhid_prior_constant)) * ((scp->tmptmp) ? scp->data.data[coor[1]].getSize() : 1);
                        printf("%e vs %e predicted\n", lls[0] - lls[2], expexp + expexp2);
                        printf("quad term: %e vs %e\n", lls[0] - lls[2] - expexp, expexp2);

                    }

                    // Step 4.1 Append LL from dropouts  TODO, ignore drop out for now
                    /*if (haschange){
                        for(i=0;i< scp->nb_maxrowh.getSize();i++{


                        }
                    }*/

                    // Step 5 select best or do not update
                    //tmp = d_scale_tmp[0] * d_scale[0] + d_scale_tmp[1] * d_scale[1];
                    //tmp += (dd_scale_tmp.data[0] * d_scale[0] * d_scale[0] + dd_scale_tmp.data[2] * d_scale[1] * d_scale[1])* 0.5 + dd_scale_tmp.data[1] * d_scale[0] * d_scale[1];
                    //printf("%e incr (vs%e)\n", lls[0] - lls[3], tmp);

                //    for(coor[0]=0;coor[0] < new_nbmax_h.getSize();coor[0]++) new_nbmax_h(coor) = scp->row_nbmax_hcol(coor);
                    if (candidates[1].getSize() == 0){
                         j = (lls[0] > lls[2]) ? 0 : 2;
                    }else{ // lls[0], lls[1] lls[3]
                        // worrying about dropouts!
                        lld.toZero();
                        acoor[1] = curcol.deref_key(0);
                        for(acoor[0] =0;acoor[0]< scp->target.par_RM.sizes[0]; acoor[0]++){coor[0] = acoor[0];
                            lld[2] += scp->col_nbmax_hrow(coor) * log(scp->target.par_D(acoor));
                            lld[2] += (scp->nb_maxrowh[coor[0]] - scp->col_nbmax_hrow(coor)) * log(1.0 - scp->target.par_D(acoor));
                        }
                        acoor[0] = candidates[0].deref_key(0);
                        for(acoor[0] =0;acoor[0]< scp->target.par_RM.sizes[0]; acoor[0]++){coor[0] = acoor[0];
                            lld[0] += scp->col_nbmax_hrow(coor) * log(scp->target.par_D(acoor));
                            lld[0] += (scp->nb_maxrowh[coor[0]] - scp->col_nbmax_hrow(coor)) * log(1.0 - scp->target.par_D(acoor));
                        }
                        if (candidates[1].deref_key(0) == curcol.deref_key(0)) lld[1] = lld[2];
                        else{
                            acoor[0] = candidates[1].deref_key(0);
                            for(acoor[0] =0;acoor[0]< scp->target.par_RM.sizes[0]; acoor[0]++){coor[0] = acoor[0];
                                lld[1] += scp->col_nbmax_hrow(coor) * log(scp->target.par_D(acoor));
                                lld[1] += (scp->nb_maxrowh[coor[0]] - scp->col_nbmax_hrow(coor)) * log(1.0 - scp->target.par_D(acoor));
                            }
                        }

                        j = (lls[0] + lld[0] > lls[2] + lld[2]) ? 0 : 2;
                        if ((lls[1] + ((lld[0] > lld[1]) ? lld[0]  : lld[1])) > lls[j] + lld[j]) {
                            if (lld[0] > lld[1]) candidates[1].swapEntries(0,1);
                            j=1;
                        }
                        if ((j != 2)&&(candidates[j].deref_key(0) != curcol.deref_key(0))){
                            // registering change...
                            nbmax_diff[candidates[j].deref_key(0)]++;
                            nbmax_diff[curcol.deref_key(0)]--;
                            for(i = 0; i < scp->data.data[coor[1]].getSize();i++){
                                acoor[0] = scp->target.rowhid.data[ scp->data.data[coor[1]].deref_key(i)].deref_key(0);
                                acoor[1] = curcol.deref_key(0);
                                nonzero_diff(acoor)--;
                                acoor[1] = candidates[j].deref_key(0);
                                nonzero_diff(acoor)++;
                            }
                        }
                    }

                    if (j != 2) break;
                    if (insist++ > 10) break;

                    if (auto ite = curcol.getIterator()) do{
                        candidates[0][ite()] += *ite;
                    }while(ite++);

                    if (auto ite = candidates[0].getIterator()) do{
                        *ite *= 0.5;
                    }while(ite++);
                    d_scale *= 0.5;



                    candidates[0].makeMaximumAsFirst();
                    if (candidates[0].deref_key(0) == curcol.deref_key(0)) candidates[1].toMemfree();

                } // insist loop end
                nb_insist += insist;

                // Step 6 save scope for parameter update
                if (j == 2) { // failed to improve...
                    //printf("col %i failed  %e %e\n", coor[1], scp->target.par_colscale[coor[1]][0], scp->target.par_colscale[coor[1]][1] );

                    candidates[0] = curcol;
                    scp->newhidden.memmoveColumn(coor[1], candidates[0]);
                    scp->newscale[coor[1]] = scp->target.par_colscale[coor[1]];
                    nb_fail++;
                }else{
                    scp->newhidden.memmoveColumn(coor[1], candidates[j]);
                    scp->newscale[coor[1]] = scp->target.par_colscale[coor[1]] + d_scale; // to change?
                }
                ll[1] += lls[j];
                candidates[1].toMemfree();

                future_factor.toMemfree();
                if (auto ite = scp->newhidden.data[coor[1]].getIterator()) do{
                    future_factor[ite()][0] = *ite; future_factor[ite()][1] = *ite;
                }while(ite++);

                dll_dMR2.toAddOuterProd(future_statistics[j], future_factor);


                uint32_t dashift = scp->target.par_RM.sizes[0] * scp->target.par_RM.sizes[1];
                if (auto ite = scp->newhidden.data[coor[1]].getIterator()) do{
                    if (auto ite2 = scp->newhidden.data[coor[1]].getIterator()) do{
                        if (ite2() > ite()) continue;
                        tmp = (*ite2) * (*ite);
                        uint32_t k;
                        for(i = 0,k=0; i< scp->target.par_RM.sizes[0];i++){
                            uint32_t l = i + ite() * scp->target.par_RM.sizes[0];

                            double* currow_A = dd_massive.data + (((l * (l +1)) >> 1) + ite2() * scp->target.par_RM.sizes[0]);
                            double* currow_B = dd_massive.data + ((((l + dashift) * (l + dashift +1)) >> 1) + ite2() * scp->target.par_RM.sizes[0]);
                            for(l=0;l<=i;l++,k++){
                                currow_A[l] += future_statistic_trix[j].data[k][0] * tmp;
                                currow_B[l] += future_statistic_trix[j].data[k][2] * tmp;
                                currow_B[l + dashift] += future_statistic_trix[j].data[k][1] * tmp;
                            }
                        }
                        //printf("%i is k  %i\n", k,  future_statistic_trix[j].getSize());
                    }while(ite2++);
                }while(ite++);
            }
        }
        lls[0] =0;
        for(coor[1]=rowrange[0];coor[1] < rowrange[1];coor[1]++){
            lls[0] += log(scp->target.rowhid_prior.mkInnerMult(scp->target.rowhid.getColumn(coor[1]))+ scp->target.rowhid_prior_constant) * ((scp->tmptmp) ? scp->trdata.data[coor[1]].getSize() : 1);
        }
        ll[0] += lls[0];
        ll[1] += lls[0];

    }else{ // time for rows!
        RM_deriv.setSize(scp->target.par_RM.sizes[0]);
        prior_deriv.setSize(scp->target.par_RM.sizes[0]);


        rever_buf.setSize(scp->target.colhid.getNBcols());
        dd_buf.setSize(scp->target.colhid.getNBcols());

        for(i=0;i<3;i++) {
            future_statistics[i].setSize(scp->target.par_RM.sizes[1]);
            future_statistic_trix[i].setSize(scp->target.par_RM.sizes[1]);
        }
        ExOp::toZero(sum);

        for(coor[1]=rowrange[0];coor[1] < rowrange[1];coor[1]++){
            if (use_stdout) pbar.update(((double)coor[1]-rowrange[0])/ (rowrange[1]-rowrange[0]));
            if (scp->trdata.data[coor[1]].getSize() > scp->target.colhid.getNBcols()) {printf("Holly Molly! %i and %i",scp->trdata.data[coor[1]].getSize(), scp->target.colhid.getNBcols());exit(1);}
            //printf("cur coor %i/%i\n", coor[1], scp->target.colhid.getNBcols());fflush(stdout);
            const SparseTuple<double>& curcol = scp->target.rowhid.getColumn(coor[1]);
            if (curcol.getSize() == 0) {printf("Queried collumn %i is empty... impossible!\n",coor[1]); exit(1);}
            if (!ExOp::isValid(scp->target.par_rowscale[coor[1]][0])){
                scp->newscale[coor[1]] = scp->target.par_rowscale[coor[1]];

                if (scp->trdata.getColSize(coor[1]) == 0) {candidates[0] = scp->target.rowhid.getColumn(coor[1]);scp->newhidden.memmoveColumn(coor[1], candidates[0]); continue;} // nothing to see...
                // column only has 1s, cannot be used to update R or M, drop out only, pure state guarrantied

                // put in stone... for now
                candidates[0] = scp->target.rowhid.getColumn(coor[1]);scp->newhidden.memmoveColumn(coor[1], candidates[0]);
                // put in stone... for now

                RM_deriv.toZero();

                for(coor[0]=0;coor[0]<scp->row_nbmax_hcol.dims[0];coor[0]++){acoor[1] = coor[0];
                    for(acoor[0] =0;acoor[0]< scp->target.par_RM.sizes[0]; acoor[0]++){
                        RM_deriv[acoor[0]][0] += scp->row_nbmax_hcol(coor) * log(scp->target.par_D(acoor));
                        RM_deriv[acoor[0]][0] += (scp->nb_maxcolh[coor[0]] - scp->row_nbmax_hcol(coor)) * log(1.0 - scp->target.par_D(acoor));
                    }
                }
                acoor[1] =0;
                for(acoor[0]=0 ;acoor[0]< scp->target.par_RM.sizes[0]; acoor[0]++){
                    if (RM_deriv[acoor[0]][0] > RM_deriv[acoor[1]][0]) acoor[1] = acoor[0];
                }
                candidates[0].toZero();
                candidates[0][acoor[1]] = 1.0;

                if (candidates[0].deref_key(0) != curcol.deref_key(0)){
                    nbmax_diff[candidates[0].deref_key(0)]++;
                    nbmax_diff[curcol.deref_key(0)]--;
                    for(i = 0; i < scp->trdata.data[coor[1]].getSize();i++){
                        acoor[1] = scp->target.colhid.data[scp->trdata.data[coor[1]].deref_key(i)].deref_key(0);
                        acoor[0] = curcol.deref_key(0);
                        nonzero_diff(acoor)--;
                        acoor[0] = candidates[0].deref_key(0);
                        nonzero_diff(acoor)++;
                    }
                }
                scp->newhidden.memmoveColumn(coor[1], candidates[0]);
            }else{

                prior_deriv = scp->target.rowhid_prior * curcol;
                // Step 1 cache line specific data
                curproj[0] = scp->target.par_RM.mkBackMult(curcol);

                //dll_dH.toZero();
                //dll_dC.toZero();

                // Step 2 find direction, find 2 vector for R and M
                RM_deriv.toZero();
                prior_llterm[0] = scp->target.rowhid_prior.mkInnerMult(curcol)+ scp->target.rowhid_prior_constant;
                lls[2] = log(prior_llterm[0]) * ((scp->tmptmp) ? scp->trdata.data[coor[1]].getSize() : 1);

                future_statistics[2].toZero();future_statistic_trix[2].toZero();
                d_scale.toZero();dd_scale.toZero();
               // printf("loopya\n"); fflush(stdout);
                for(i=0;i<scp->trdata.data[coor[1]].getSize();i++){
                    coor[0]= scp->trdata.data[coor[1]].deref_key(i);
                    if (!ExOp::isValid(scp->target.par_colscale[coor[0]][0])) continue; // can do better...
                    coeff = curproj[0].mkInnerProd(scp->target.colhid.getColumn(coor[0]) );
                    coeff += scp->target.par_colscale[coor[0]];
                    coeff += scp->target.par_rowscale[coor[1]];
                    deflatedNBterm_wr_f_d2(ld.evilterms, scp->trdata.data[coor[1]].deref(i), coeff[0], coeff[1]);
                    dd_scale.data[0] += ld.evilterms[3]; dd_scale.data[1] += ld.evilterms[5]; dd_scale.data[2] += ld.evilterms[4];
                    d_scale[0] += ld.evilterms[1]; d_scale[1] += ld.evilterms[2];
                    lls[2] += ld.evilterms[0];

                    rever = scp->target.par_RM * scp->target.colhid.getColumn(coor[0]);
                    if (auto ite = rever.getIterator()) do{
                        RM_deriv[ite()][0] += (*ite)[0] * ld.evilterms[1];
                        RM_deriv[ite()][1] += (*ite)[1] * ld.evilterms[2];
                    }while(ite++);


                    if (auto ite = scp->target.colhid.getColumn(coor[0]).getIterator()) do{
                        future_statistics[2][ite()][0] += ld.evilterm_alt.d[0] * (*ite);
                        future_statistics[2][ite()][1] += ld.evilterm_alt.d[1] * (*ite);
                        trix_input[0] = ld.evilterms[3] * (*ite);
                        trix_input[1] = ld.evilterms[4] * (*ite);
                        trix_input[2] = ld.evilterms[5] * (*ite) * 0.5;
                        if (auto ite2 = scp->target.colhid.getColumn(coor[0]).getIterator()) do{
                            if (ite2() > ite()) continue;
                            future_statistic_trix[2].data[ite2() + ((ite() * (ite() + 1)) >> 1)] += trix_input * (*ite2);
                        }while(ite2++);

                    }while(ite++);

                    rever_buf[i].toMemmove(rever);
                    dd_buf[i][0] = ld.evilterms[0];
                    dd_buf[i][1] = ld.evilterms[1];
                    dd_buf[i][2] = ld.evilterms[2];
                    dd_buf[i][3] = ld.evilterms[3];
                    dd_buf[i][4] = ld.evilterms[4];
                    dd_buf[i][5] = ld.evilterms[5];
                }
                //dd_scale.data[1] *= 0.5;

                ll[0] += lls[2];
                // Step 3 (line search) define set of candidates


                sum[0] = 0.0;
                prior_llterm[2] = 0.0;
                candidates[0].toZero();
                for(i=0;i< RM_deriv.getSize();i++){
                    prior_llterm[1] = RM_deriv[i][0] + RM_deriv[i][1] + prior_deriv[i] * (2.0 * ((scp->tmptmp) ? scp->trdata.data[coor[1]].getSize() : 1) / prior_llterm[0]);
                    if ((curcol.find(i) != 0xFFFFFFFF)||(prior_llterm[1] > 0.0f)) {
                        candidates[0][i] = prior_llterm[1];
                        sum[0] +=  prior_llterm[1];
                    }
                }
                if (candidates[0].getSize() < 2u) { // stuck in a corner,... that's good but...
                    candidates[0] = curcol;
                    dd_scale_tmp = dd_scale; d_scale_tmp = d_scale;
                    d_scale = dd_scale.toEigenTransform(funabsoft) * d_scale;
                    d_scale *= scp->hiddenalpha;
                }else{
                    // project onto plane
                    coor[0] =0;
                    can_backup = candidates[0];
                    do{
                        torem.toMemfree();
                        sum[0] /= candidates[0].getSize();
                        if (auto ite = candidates[0].getIterator()) do{
                            if (((*ite) > sum[0])||(curcol.find(ite()) != 0xFFFFFFFF)) (*ite) -= sum[0];
                            else torem.push_back(ite());
                        }while(ite++);
                        if (torem.getSize() == 0u) {
                            sum[0] = 0.0;
                            if (auto ite = candidates[0].getIterator()) do{sum[0] += (*ite);}while(ite++);
                            if (fabs(sum[0]) < 0.0001) break;
                            if ((coor[0]++) == 5) {
                                printf("warning: got stuck... catastrophy cancellations\n");
                                can_backup.show();
                                curcol.show();
                                coor[0] = 0;
                                if (auto ite = candidates[0].getIterator()) do{if ((*ite)>= 0) coor[0]++;}while(ite++);
                                if (auto ite = candidates[0].getIterator()) do{(*ite) = ((*ite)>= 0) ?  candidates[0].getSize() - coor[0]: coor[0];}while(ite++);
                                if (auto ite = candidates[0].getIterator()) do{(*ite) /=  candidates[0].getSize();}while(ite++);
                                exit(1);
                                break;
                            }
                        }else{

                            if (candidates[0].getSize() - torem.getSize() <= 1u){ // would leave nothing... try the highest second?
                                if (candidates[0].getSize() == torem.getSize()) {printf("Impossibleble!\n"); exit(1);}/*{ // wow, rounding error causes this
                                    acoor[0] = 0;
                                    for(acoor[1]=1;acoor[1] < torem.getSize();acoor[1]++) if (candidates[0][torem[acoor[0]]] < candidates[0][torem[acoor[1]]]) acoor[0] = acoor[1];
                                    torem.pop_swap(acoor[0]);
                                }*/

                                acoor[0] = 0;
                                for(acoor[1]=1;acoor[1] < torem.getSize();acoor[1]++) if (candidates[0][torem[acoor[0]]] < candidates[0][torem[acoor[1]]]) acoor[0] = acoor[1];
                                for(i=0;i<torem.getSize();i++) {
                                    if (i == acoor[0]) continue;
                                    j = candidates[0].find(torem[i]);
                                    if (j == 0xFFFFFFFF) {printf("erasing unexisting... !\n"); exit(1);}
                                    candidates[0].erase_from_iterator(j);
                                }
                                if (candidates[0].getSize() != 2) {printf("2 should remain. got %i\n", candidates[0].getSize()); exit(1);}
                                if (candidates[0].deref_key(0) == torem[acoor[0]]){
                                    candidates[0].deref(0) = 1.0;
                                    candidates[0].deref(1) = -1.0;
                                }else{
                                    candidates[0].deref(0) = -1.0;
                                    candidates[0].deref(1) = 1.0;
                                }
                                break;
                            }
                            for(i=0;i<torem.getSize();i++) {
                                j = candidates[0].find(torem[i]);
                                if (j == 0xFFFFFFFF) {printf("trying to remove un-existing!\n"); exit(1);}
                                candidates[0].erase_from_iterator(j);
                            }
                            sum[0] =0.0;
                            if (auto ite = candidates[0].getIterator()) do{sum[0] += *ite;}while(ite++);
                        }
                    }while(true);


                    if (auto ite = candidates[0].getIterator()){
                        while((*ite) > 0.0) ite++;
                        sum[0] = curcol[dir = ite()] / (*ite);

                        while(ite++){
                            if (*ite > 0.0) continue;
                            sum[1] = curcol[ite()] / (*ite);
                            if (sum[1] > sum[0]) {sum[0] = sum[1]; dir = ite();}
                        }
                        sum[0] = -sum[0]; // w has to be in 0-sum[0] range

                        tosolve_num2[0] = (scp->target.rowhid_prior.mkInnerMult(candidates[0],curcol)+ scp->target.rowhid_prior_constant) * ( 2.0 * ((scp->tmptmp) ? scp->trdata.data[coor[1]].getSize() : 1)  / prior_llterm[0]);
                        tosolve_num2[1] = d_scale[0];
                        tosolve_num2[2] = d_scale[1];
                        tosolve_den2.data[1] = 0.0; tosolve_den2.data[3] = 0.0;
                        tosolve_den2.data[2] = dd_scale.data[0];
                        tosolve_den2.data[5] = dd_scale.data[2];
                        tosolve_den2.data[4] = dd_scale.data[1];
                        tosolve_den2.data[0] = (scp->target.rowhid_prior.mkInnerMult(candidates[0])+ scp->target.rowhid_prior_constant) * ( 2.0 * ((scp->tmptmp) ? scp->trdata.data[coor[1]].getSize() : 1) / prior_llterm[0]);

                        for(i=0;i<scp->trdata.data[coor[1]].getSize();i++){
                            coor[0]= scp->trdata.data[coor[1]].deref_key(i);
                            if (!ExOp::isValid(scp->target.par_colscale[coor[0]])) continue;
                            coeff = rever_buf[i].mkInnerProd(candidates[0]);
                            tosolve_num2[0] += coeff[0] * dd_buf[i][1] + coeff[1] * dd_buf[i][2];
                            tosolve_den2.data[0] += dd_buf[i][3] * coeff[0] * coeff[0] + dd_buf[i][4] * coeff[1] * coeff[1] + dd_buf[i][5] * coeff[0] * coeff[1];
                            tosolve_den2.data[1] += dd_buf[i][3] * coeff[0] + dd_buf[i][5] * coeff[1] * 0.5;
                            tosolve_den2.data[3] += dd_buf[i][4] * coeff[1] + dd_buf[i][5] * coeff[0] * 0.5;
                        }

                        tosolve_den2_tmp = tosolve_den2;
                        tosolve_num2_tmp = tosolve_num2;
                        tosolve_den2.toEigenTransform(funabsoft);
                        tosolve_num2 = tosolve_den2 * tosolve_num2;
                        tosolve_num2 *= scp->hiddenalpha;

                        expexp = tosolve_num2[0] * tosolve_num2_tmp[0] + tosolve_num2[1] * tosolve_num2_tmp[1] + tosolve_num2[2] * tosolve_num2_tmp[2];
                        expexp2 = tosolve_den2_tmp.data[0] * tosolve_num2[0] * tosolve_num2[0] * 0.5;
                        expexp2 += tosolve_den2_tmp.data[2] * tosolve_num2[1] * tosolve_num2[1] * 0.5;
                        expexp2 += tosolve_den2_tmp.data[5] * tosolve_num2[2] * tosolve_num2[2] * 0.5;
                        expexp2 += tosolve_den2_tmp.data[1] * tosolve_num2[0] * tosolve_num2[1];
                        expexp2 += tosolve_den2_tmp.data[3] * tosolve_num2[0] * tosolve_num2[2];
                        expexp2 += tosolve_den2_tmp.data[4] * tosolve_num2[1] * tosolve_num2[2];

                        if (tosolve_num2[0] <= 0.0){ // went the other way... just update scale
                            if (tosolve_num2[0] == 0.0){
                                candidates[0] = curcol;
                            }else{
                                candidates[0] = curcol;
                                dd_scale_tmp = dd_scale; d_scale_tmp = d_scale;
                                d_scale = dd_scale.toEigenTransform(funabsoft) * d_scale;
                                d_scale *= scp->hiddenalpha;
                            }
                        }else if (tosolve_num2[0] >= sum[0]){ // too far!, find solution on boundary...
                            d_scale[0] = tosolve_num2_tmp[1] + sum[0] * tosolve_den2.data[1];
                            d_scale[1] = tosolve_num2_tmp[2] + sum[0] * tosolve_den2.data[3];
                            d_scale = dd_scale.toEigenTransform(funabsoft) * d_scale;
                            d_scale *= scp->hiddenalpha;
                            nb_funky++;
                            candidates[1] = candidates[0];
                            if (auto ite = candidates[0].getIterator()) do{
                                (*ite) = (*ite) * sum[0]  + (((j = curcol.find(ite())) == 0xFFFFFFFF) ? 0.0: curcol.deref(j));
                            }while(ite++);

                            if (auto ite = curcol.getIterator()) do{
                                if (candidates[0].find(ite())== 0xFFFFFFFF) candidates[0][ite()] = *ite;
                            }while(ite++);

                            j = candidates[0].find(dir);
                            if (j == 0xFFFFFFFF) {printf("trying to remove un-existing!\n");exit(1);}
                            candidates[0].erase_from_iterator(j);
                            sum[1] =0.0;
                            if (auto ite = candidates[0].getIterator()) {
                                do{sum[1] += *ite;}while(ite++);
                                if (fabs(sum[1] - 1.0) > 0.1) {printf("bound proj dir %i c... hidden state norm is %e! %e is fact\n", dir, sum[1], sum[0]);
                                candidates[0].show();
                                candidates[1].show();
                                curcol.show();
                                LFH_ALIVE;exit(1);}
                            }else {printf("got empty candidate!\n"); exit(1);}
                            candidates[1].toMemfree();
                        }else{ // all set!
                            d_scale[0] = tosolve_num2[1];
                            d_scale[1] = tosolve_num2[2];
                            candidates[1] = candidates[0];
                            if (auto ite = candidates[0].getIterator()) do{
                                (*ite) = (*ite) * tosolve_num2[0] + (((j = curcol.find(ite())) == 0xFFFFFFFF) ? 0.0: curcol.deref(j));
                            }while(ite++);
                            if (auto ite = curcol.getIterator()) do{
                                if (candidates[0].find(ite())== 0xFFFFFFFF) candidates[0][ite()] = *ite;
                            }while(ite++);
                            sum[1]= 0.0f;
                            if (auto ite = candidates[0].getIterator()) {
                                do{sum[1] += *ite;}while(ite++);
                                if (fabs(sum[1] - 1.0) > 0.1) {printf("normal proj rowr... hidden state norm is %e! %e is fact\n", sum[1], tosolve_num2[0]);
                                candidates[1].show();
                                candidates[0].show();
                                curcol.show();
                                exit(1);}
                            }else {printf("got empty candidate!\n"); exit(1);}
                            candidates[1].toMemfree();
                        }
                        candidates[0].makeMaximumAsFirst();
                        if ((candidates[0].deref_key(0) != curcol.deref_key(0))&&(!scp->ignore_drop)){
                            // major change, add *last* boundary on the search space!
                            tosolve_num2[2] = ((i = curcol.find(candidates[0].deref_key(0))) == 0xFFFFFFFF) ? 0.0 : curcol.deref(i);
                            if (auto ite = curcol.getIterator()) {
                                tosolve_num2[0] = (*ite) - tosolve_num2[2]; // guarranied to be positive
                                tosolve_num2[1] = candidates[0].deref(0) -  (((i = candidates[0].find(ite())) == 0xFFFFFFFF) ? 0.0 : candidates[0].deref(i));
                                sum[1] = tosolve_num2[1] / (tosolve_num2[0] + tosolve_num2[1]);
                                j = 0;
                                while(ite++){
                                    if ((tosolve_num2[0] =  (*ite) - tosolve_num2[2]) <= 0) continue;
                                    tosolve_num2[1] = candidates[0].deref(0) -  (((i = candidates[0].find(ite())) == 0xFFFFFFFF) ? 0.0 : candidates[0].deref(i));
                                    sum[0] = tosolve_num2[1] / (tosolve_num2[0] + tosolve_num2[1]);
                                    if (sum[0] < sum[1]) {j = ite.getOffset(); sum[1] = sum[0];}
                                }
                                sum[0] = candidates[0].deref(0) * (1.0 - sum[1]) + sum[1] * (((i = curcol.find(candidates[0].deref_key(0))) == 0xFFFFFFFF) ? 0.0 : curcol.deref(i));
                                candidates[1][curcol.deref_key(j)] = sum[0];
                                candidates[1][candidates[0].deref_key(0)] = sum[0];
                                if (auto ite = curcol.getIterator()) do {
                                    if (candidates[1].find(ite()) == 0xFFFFFFFF) candidates[1][ite()] = (*ite) * sum[1] + (1.0 - sum[1]) * (((i = candidates[0].find(ite())) == 0xFFFFFFFF) ? 0.0 : candidates[0].deref(i));
                                }while(ite++);
                                if (auto ite = candidates[0].getIterator()) do {
                                    if (candidates[1].find(ite()) == 0xFFFFFFFF) candidates[1][ite()] = (*ite) * (1.0 - sum[1]);
                                }while(ite++);
                            }

                            if (auto ite = candidates[1].getIterator()) {
                                sum[0] =0.0;
                                do{
                                    sum[0] += *ite;
                                }while(ite++);
                                if (fabs(sum[0] - 1.0) > 0.1){printf("hidded state norm cancan1 is %e!\n", sum[0]);
                                    candidates[1].show();
                                    printf("frac %e\n", sum[1]);
                                    candidates[0].show();
                                    curcol.show();
                                exit(1);}
                            }else {printf("got empty candidate!\n"); exit(1);}

                        }
                    }else{
                        // stuck in corner again...
                        candidates[0] = curcol;
                        dd_scale_tmp = dd_scale; d_scale_tmp = d_scale;
                        d_scale = dd_scale.toEigenTransform(funabsoft) * d_scale;
                        d_scale *= scp->hiddenalpha;
                    }
                }
                sum[0] = 0.0f;
                if (auto ite = candidates[0].getIterator()) {
                    do{
                        sum[0] += *ite;
                    }while(ite++);
                    if (fabs(sum[0] - 1.0) > 0.1){printf("hidded state norm is %e!\n", sum[0]);
                        candidates[0].show();
                        candidates[1].show();
                    exit(1);}
                }else {printf("got empty candidate!\n"); exit(1);}
                double tmp;

                //printf("cur coor %i/%i\n", coor[1], scp->target.colhid.getNBcols());fflush(stdout);

                coeff = scp->target.par_rowscale[coor[1]];
                coeff += d_scale;
                if ((tmp = (coeff[0]*coeff[0] + coeff[1]*coeff[1])) > 100.0){ // too far, send to boundary
                    tmp = pow(tmp, -0.5);
                    d_scale[0] = coeff[0] * tmp;
                    d_scale[1] = coeff[1] * tmp;
                    d_scale -= scp->target.par_rowscale[coor[1]];
                }

                insist =0;
                while(true){

                    // Step 4 Evaluate Quantidates
                    for(j=0;j<2;j++){
                        future_statistics[j].toZero();future_statistic_trix[j].toZero();
                        curproj[j] = scp->target.par_RM.mkBackMult(candidates[j]);
                        lls[j] = log(scp->target.rowhid_prior.mkInnerMult(candidates[j])+ scp->target.rowhid_prior_constant) * ((scp->tmptmp) ? scp->trdata.data[coor[1]].getSize() : 1) ;
                        if ((candidates[1].getSize() == 0)||(insist > 0)) break;
                    }
                                 //   printf("loopyb\n"); fflush(stdout);
                    for(i=0;i<scp->trdata.data[coor[1]].getSize();i++){
                        coor[0]= scp->trdata.data[coor[1]].deref_key(i);
                        if (!ExOp::isValid(scp->target.par_colscale[coor[0]])) continue; // can do better...
                        for(j=0;j<2;j++){
                            coeff = curproj[j].mkInnerProd(scp->target.colhid.getColumn(coor[0]));
                            coeff += scp->target.par_colscale[coor[0]];
                            coeff += scp->target.par_rowscale[coor[1]];
                            coeff += d_scale;
                            deflatedNBterm_wr_f_d2(ld.evilterms, scp->trdata.data[coor[1]].deref(i), coeff[0], coeff[1]);
                           if (auto ite = scp->target.colhid.getColumn(coor[0]).getIterator()) do{
                                future_statistics[j][ite()][0] += ld.evilterm_alt.d[0] * (*ite);
                                future_statistics[j][ite()][1] += ld.evilterm_alt.d[1] * (*ite);

                                trix_input[0] = ld.evilterms[3] * (*ite);
                                trix_input[1] = ld.evilterms[4] * (*ite);
                                trix_input[2] = ld.evilterms[5] * (*ite) * 0.5;
                                if (auto ite2 = scp->target.colhid.getColumn(coor[0]).getIterator()) do{
                                    if (ite2() > ite()) continue;
                                    future_statistic_trix[j].data[ite2() + ((ite() * (ite() + 1)) >> 1)] += trix_input * (*ite2);
                                }while(ite2++);

                            }while(ite++);

                            lls[j] += ld.evilterms[0];
                            if ((candidates[1].getSize() == 0)||(insist > 0)) break;
                        }
                    }
                    //      printf("loopyb done\n"); fflush(stdout);

                    /*if (!checksafe) nb_funky++;
                    else {ld.evilterms[1]= (1.0 - (lls[0] - lls[3]) / prediction); accuracy += ld.evilterms[1]*ld.evilterms[1];}*/


                    // Step 4.1 Append LL from dropouts  TODO, ignore drop out for now
                    /*if (haschange){
                        for(i=0;i< scp->nb_maxrowh.getSize();i++{


                        }
                    }*/

                    // Step 5 select best or do not update
                   // tmp = d_scale_tmp[0] * d_scale[0] + d_scale_tmp[1] * d_scale[1];
                   // tmp += (dd_scale_tmp.data[0] * d_scale[0] * d_scale[0] + dd_scale_tmp.data[2] * d_scale[1] * d_scale[1])* 0.5 + dd_scale_tmp.data[1] * d_scale[0] * d_scale[1];
                   // printf("%e incr (vs%e)\n", lls[0] - lls[3], tmp);

                    if (candidates[1].getSize() == 0){
                         j = (lls[0] > lls[2]) ? 0 : 2;
                    }else{
                        // worrying about dropouts!
                        lld.toZero();
                        acoor[0] = curcol.deref_key(0);
                        for(acoor[1] =0;acoor[1]< scp->target.par_RM.sizes[1]; acoor[1]++){coor[0] = acoor[1];
                            lld[2] += scp->row_nbmax_hcol(coor) * log(scp->target.par_D(acoor));
                            lld[2] += (scp->nb_maxcolh[coor[0]] - scp->row_nbmax_hcol(coor)) * log(1.0 - scp->target.par_D(acoor));
                        }
                        acoor[0] = candidates[0].deref_key(0);
                        for(acoor[1] =0;acoor[1]< scp->target.par_RM.sizes[1]; acoor[1]++){coor[0] = acoor[1];
                            lld[0] += scp->row_nbmax_hcol(coor) * log(scp->target.par_D(acoor));
                            lld[0] += (scp->nb_maxcolh[coor[0]] - scp->row_nbmax_hcol(coor)) * log(1.0 - scp->target.par_D(acoor));
                        }
                        if (candidates[1].deref_key(0) == curcol.deref_key(0)) lld[1] = lld[2];
                        else{
                            acoor[0] = candidates[1].deref_key(0);
                            for(acoor[1] =0;acoor[1]< scp->target.par_RM.sizes[1]; acoor[1]++){coor[0] = acoor[1];
                                lld[1] += scp->row_nbmax_hcol(coor) * log(scp->target.par_D(acoor));
                                lld[1] += (scp->nb_maxcolh[coor[0]] - scp->row_nbmax_hcol(coor)) * log(1.0 - scp->target.par_D(acoor));
                            }
                        }
                        j = (lls[0] + lld[0] > lls[2] + lld[2]) ? 0 : 2;
                        if ((lls[1] + ((lld[0] > lld[1]) ? lld[0]  : lld[1])) > lls[j] + lld[j]){
                            if (lld[0] > lld[1]) candidates[1].swapEntries(0,1);
                            j=1;
                        }
                        if ((j != 2)&&(candidates[j].deref_key(0) != curcol.deref_key(0))){
                            // registering change...
                            nbmax_diff[candidates[j].deref_key(0)]++;
                            nbmax_diff[curcol.deref_key(0)]--;
                            for(i = 0; i < scp->trdata.data[coor[1]].getSize();i++){
                                acoor[1] = scp->target.colhid.data[scp->trdata.data[coor[1]].deref_key(i)].deref_key(0);
                                acoor[0] = curcol.deref_key(0);
                                nonzero_diff(acoor)--;
                                acoor[0] = candidates[j].deref_key(0);
                                nonzero_diff(acoor)++;
                            }
                        }
                    }

                    if (j != 2) break;
                    if (insist++ > 10) break;

                    if (auto ite = curcol.getIterator()) do{
                        candidates[0][ite()] += *ite;
                    }while(ite++);

                    if (auto ite = candidates[0].getIterator()) do{
                        *ite *= 0.5;
                    }while(ite++);
                    d_scale *= 0.5;


                    candidates[0].makeMaximumAsFirst();
                    if (candidates[0].deref_key(0) == curcol.deref_key(0)) candidates[1].toMemfree();

                }
                nb_insist += insist;
        //printf("curss coor %i/%i\n", coor[1], scp->target.colhid.getNBcols());fflush(stdout);

                // Step 6 save scope for parameter update
                if (j == 2) { // failed to improve...
                    //printf("row %i failed %e %e\n", coor[1], scp->target.par_rowscale[coor[1]][0], scp->target.par_rowscale[coor[1]][1] );
                    candidates[0] = curcol;
                    scp->newhidden.memmoveColumn(coor[1], candidates[0]);
                    scp->newscale[coor[1]] = scp->target.par_rowscale[coor[1]];
                    nb_fail++;
                }else{
                    scp->newhidden.memmoveColumn(coor[1], candidates[j]);
                    scp->newscale[coor[1]] = scp->target.par_rowscale[coor[1]] + d_scale;
                }
                ll[1] += lls[j];
                candidates[1].toMemfree();

                future_factor.toMemfree();
                if (auto ite = scp->newhidden.data[coor[1]].getIterator()) do{
                    future_factor[ite()][0] = *ite; future_factor[ite()][1] = *ite;
                }while(ite++);
                dll_dMR2.toAddBackOuterProd(future_statistics[j], future_factor);


                uint32_t dashift = scp->target.par_RM.sizes[0] * scp->target.par_RM.sizes[1];
                uint32_t k;
                for(i = 0,k=0; i< scp->target.par_RM.sizes[1];i++){
                    for(uint32_t l=0;l<=i;l++,k++){
                        if (auto ite = scp->newhidden.data[coor[1]].getIterator()) do{
                            uint32_t m = ite() + i * scp->target.par_RM.sizes[0];
                            double* currow_A = dd_massive.data + (((m * (m +1)) >> 1) + l * scp->target.par_RM.sizes[0]);
                            double* currow_B = dd_massive.data + ((((m + dashift) * (m + dashift +1)) >> 1) + l * scp->target.par_RM.sizes[0]);
                            if (auto ite2 = scp->newhidden.data[coor[1]].getIterator()) do{
                                if (ite2() > ite()) continue;
                                double tmp = (*ite) * (*ite2);
                                currow_A[ite2()] += future_statistic_trix[j].data[k][0] * tmp;
                                currow_B[ite2()] += future_statistic_trix[j].data[k][2] * tmp;
                                currow_B[ite2() + dashift] += future_statistic_trix[j].data[k][1] * tmp;
                            }while(ite2++);
                        }while(ite++);
                    }
                }
            }
        }
        lls[0] =0;
        for(coor[1]=colrange[0];coor[1] < colrange[1];coor[1]++){
            lls[0] += log(scp->target.colhid_prior.mkInnerMult(scp->target.colhid.getColumn(coor[1]))+ scp->target.colhid_prior_constant) * ((scp->tmptmp) ? scp->data.data[coor[1]].getSize() : 1);
        }
        ll[0] += lls[0];
        ll[1] += lls[0];
    }
    if (use_stdout) {
        pbar.finish();
        printf("%e accuracy %c %e \n", alpha, (scp->istimefor_columns) ? 'c': 'r', accuracy / ((scp->istimefor_columns) ? colrange[1] - nb_funky : rowrange[1] - nb_funky));
    }
    return 1000;
}

#undef LFHTEMP




template<Tuple_flag TF> void MatrixDistribution::expect_routine(Tuple<double, 0u, TF> &_out, const Tuple<double, 0u, TF> &_in, const double* _stat_vector, const double difsum, const Trianglix<double> innersum, const Trianglix<double> &mat_U, const uint32_t _nb, Tuple<double> &der, Tuple<double> &dir, Trianglix<double> &dder){
    // this step needs to be guarantied to increase the likelihood! And it bounces around the hyper-tetrahedron boundary.
    // f(Z) = (Z-U_p) M_p (Z-U_p) + ...
    // f'(Z) = 2M_p (Z-U_p) +  (_nb U Z + _stat_vector) / ZtUZ + difsum / ZtUZ^2
    // f"(Z) = 2M_p + ...
    double inner = mat_U.Xformed_inner_product(_in);
    unsigned int i;
    // Step 1: get direction for Newton step
    for(i=0;i<der.getSize();i++){
        der += mat_U * _in;

    }

    dir = dder.mkInvDivi(der);

    // Step 2: align direction within domain boundaries

    double sum = 1.0f;
    double dotprod = 0.0f;
    unsigned int nbvalid =der.getSize();
    for(i=0;i<der.getSize();i++){

        if (_in[i] == 0.0f){
            nbvalid--;
            if (dir[i] < 0.0f) dir[i] =0.0f; // do not go in negative...
        }else {sum -= _in[i];}
        dotprod += dir[i];
    }
    if ((sum <= 2.2204460492503131e-16 * der.getSize())&&(dotprod > 0.0f)) { // make the projection in orthogonal direction 0
        for(i=0;i<der.getSize();i++){
            if (_in[i] != 0.0f) dir[i] -= dotprod / nbvalid;
        }
    }

    // Step 3: line search,

    //dir = ddir * dir;

    double proj =0;


    for(i=0;i<der.getSize();i++) proj +=  dir[i];
    proj = (proj > 0) ?  proj : 0;
    int best = der.getSize();
    for(i=0;i<der.getSize();i++) {
        if ((dir[i] < 0)&&( proj < -dir[i])) {best = i; proj = -dir[i];}
    }




}
template<class C, unsigned int DIMS> DataGrid<double,2> MetaPixels::performEM(unsigned int nbstep, const DataGrid<C, DIMS> &obs, ClassifierV<C> &model, const DataGrid<unsigned int, DIMS> &meta, unsigned int nbmetapix, unsigned int first_val){
	Vector< C > *pixlists = new Vector< C >[nbmetapix];
	unsigned int i,j;

	typename DataGrid<C, DIMS>::KeyIterator ite = obs.getKeyIterator();

	if (ite.first()) do{
		j = meta(ite());
		if ((j >= first_val)&&(j < first_val + nbmetapix)) pixlists[j -  first_val].push_back(obs(ite()));
	} while ( ite.next());


	unsigned int nbstates = model.nbstates();


	double *likeli = new double[nbstates];
	double* prior_freq = new double[nbstates];
	double* prior_freq_learn = new double[nbstates];

	for(j=0;j< nbstates; j++) prior_freq[j] = 1.0f / nbstates;



	int r;
	unsigned int x,y;

	unsigned int totnbpix =0;for(y = 0;y< nbmetapix;y++) totnbpix += pixlists[y].size();

	double LLhood;
	double LLhood_old;

	double alpha = 1.0f; //metapix
	double sum,saf;


	bool rand_init = false;
	for(i=0;i<nbstep;i++){

		model.EMinit();
		memset(likeli,'\0',sizeof(double)*nbstates);
		memset(prior_freq_learn,'\0',sizeof(double)*nbstates);
		if ((i==0)&&(rand_init)){
			for(y = 0;y< nbmetapix;y++) {
				r = rand() % nbstates;
				for(j=0;j< nbstates; j++) likeli[j] = (r == j) ? 1.0f : 0.0f;
			for(x = 0;x< pixlists[y].size();x++) model.EMregist(pixlists[y][x] , likeli);
			}
		}else{
			for(y = 0;y< nbmetapix;y++) {
				memset(likeli,'\0',sizeof(double)*nbstates);
				for(x = 0;x< pixlists[y].size();x++) for(j=0;j< nbstates; j++) likeli[j] += model[j]->LL(pixlists[y][x]);
				saf = likeli[0];
				for(j=1;j< nbstates; j++) if (saf < likeli[j]) saf = likeli[j];
				sum =0;
				for(j=0;j< nbstates; j++) {likeli[j] = prior_freq[j] * exp(likeli[j] - saf); sum += likeli[j];}
				for(j=0;j< nbstates; j++) likeli[j] /= sum;
				for(x = 0;x< pixlists[y].size();x++) model.EMregist(pixlists[y][x] , likeli);
				for(j=0;j< nbstates; j++) prior_freq_learn[j] += likeli[j] * pixlists[y].size();
			}
		}
		LLhood = model.EMfinit();
		for(j=0;j< nbstates; j++) prior_freq[j] = (ExCo<double>::epsilon() + prior_freq_learn[j]) / ( (ExCo<double>::epsilon() *nbstates) + totnbpix);

		if (((LLhood_old > LLhood)||(!ExCo<double>::isValid(LLhood)))&&(i != 0)) {
			printf("Failed: Log-likelihood = %e\n", LLhood);
			alpha *= 0.1f;
			model.EMaccept(false);
		}else{
			printf("Step %i: Log-likelihood = %e\n", i , LLhood);
			LLhood_old = LLhood;
			alpha = 1.0f;
			model.EMaccept();

		}
		model.EMnext(alpha);
	}
	model.EMclear();

	DataGrid<double,2> fout;
	unsigned int coor[2];
	coor[0] = nbstates;
	coor[1] = nbmetapix;
	fout.setSizes(coor);


	for(coor[1] = 0;coor[1]< nbmetapix;coor[1]++) {
		memset(likeli,'\0',sizeof(double)*nbstates);
		for(x = 0;x< pixlists[coor[1] ].size();x++) for(j=0;j< nbstates; j++) likeli[j] += model[j]->LL(pixlists[coor[1]][x]);
		saf = likeli[0];
		for(j=1;j< nbstates; j++) if (saf < likeli[j]) saf = likeli[j];
		sum =0;
		for(j=0;j< nbstates; j++) {likeli[j] = prior_freq[j] * exp(likeli[j] - saf); sum += likeli[j];}
		for(coor[0]=0;coor[0]< nbstates; coor[0]++) fout(coor) = likeli[ coor[0] ] / sum;
	}

	delete[](likeli);
	delete[](prior_freq);
	delete[](prior_freq_learn);

	delete[](pixlists);
	return fout;
}



#undef LFHTEMP
#define LFHTEMP template<class C, unsigned int SIZE>


LFHTEMP typename ExCo<C>::LMUL_TYPE GaussianObj<C,SIZE>::computeProj(const Tuple<C,SIZE> &) const{

}

LFHTEMP void GaussianObj<C,SIZE>::operator()(double &fout , Tuple<C,SIZE> &fin) const{
    typename ExCo<C>::LMUL_TYPE mag = computeProj(fin);
    fout = factor * exp(-ExOp::norm(mag));
}

LFHTEMP double GaussianObj<C,SIZE>::LL(const C& fin) const{
    typename ExCo<C>::LMUL_TYPE mag = computeProj(fin);
    return log(factor) -ExOp::norm(mag);
}

LFHTEMP void GaussianObj<C,SIZE>::EMinit(){
 //   delete[](dascope); dascope = new GOscope();
}
LFHTEMP void GaussianObj<C,SIZE>::EMAlphainit(double a){
 //   delete[](dascope); dascope = new GOscope();
}

LFHTEMP double GaussianObj<C,SIZE>::EMfinit(){
	mean = scope.getMean();
	ihvar = scope.cov.inverse() * -0.5f;
    return(0.0f);
}

LFHTEMP void GaussianObj<C,SIZE>::EMregist(const Tuple<C,SIZE> &instance, const double prob){
	scope += GaussElem<C,SIZE>(instance, prob);
}

LFHTEMP GaussianObj<C,SIZE>* GaussianObj<C,SIZE>::EMnext(double alpha) const{
    return new GaussianObj<C,SIZE>(*this);
}

LFHTEMP void GaussianObj<C,SIZE>::EMclear(){
//delete[](dascope);dascope = NULL;
}

LFHTEMP void GaussianObj<C,SIZE>::show(FILE *f, unsigned int level)const{
 //   fprintf(f,"Gaussian Object: mean = "); ExOp::show(mean, f, level+1);
 //   fprintf(f,"\tProject = ");ExOp::show(proj, f, level+1);
 //   fprintf(f,"\n");
}


#undef LFHTEMP
#define LFHTEMP template<unsigned int NBPARAM, unsigned int NBINPUT>


LFHTEMP	void GaussianProcessProbSimple<NBPARAM,NBINPUT>::EMinit(){

	noise.EMinit();
	EM_obs.toMemfree();

/*
	if (resample_to_fixgrid){

		Tuple<double, NBINPUT> daz; ExOp::toZero(daz);
		for(unsigned int i=0;i<registered.size();i++){registered[i].first.d = daz; registered[i].second = 0.01f;}
	}else{
		registered.clear();
	}
*/

}

LFHTEMP	void GaussianProcessProbSimple<NBPARAM,NBINPUT>::EMAlphainit(double a){
	noise.EMAlphainit(a);
    noise.EMinit();
	EM_obs.clear();
}


LFHTEMP	void GaussianProcessProbSimple<NBPARAM,NBINPUT>::EMregist(const KeyElem< Tuple<double, NBPARAM> , Tuple<double, NBINPUT> > &data,  const double prob){
    if (resample_to_fixgrid){
    }else EM_obs.push_back(KeyElem< double, KeyElem< Tuple<double, NBPARAM> , Tuple<double, NBINPUT> > >(prob, data));

	noise.EMregist(data.d, prob);
}


LFHTEMP	void GaussianProcessProbSimple<NBPARAM,NBINPUT>::EMclear(){
	EM_obs.clear();
}


LFHTEMP double GaussianProcessProbSimple<NBPARAM,NBINPUT>::setDerivative_routine(double* i_deriv, bool showmat){
    Trianglix<double, 0u> damatrix;
    Trianglix<double, 0u> damatrix_nonoise;

	Trianglix<double, 0u> daderivative;
	Trianglix<double, 0u> dainverse;
    Trianglix<double, 0u> dainner;

	inner_param.setSize(EM_obs.size());

	Tuple<double,0u> datuple;
Matrix<double> test;


    double ttmptmp,tmp_dev;

 //   DataGrid<double, 2> da_vectors;
   // DataGrid<double, 1> inner_prods;
 //   unsigned int coor[2];
 //   coor[0] = 1;
    unsigned int s = EM_obs.size();
 //   coor[1] = s;
 //   da_vectors.setSizes(coor);
    unsigned int i,j,k;
 //   printf("%e\t%e\t%e\n",scope[0],scope[1],scope[2]);
				damatrix_nonoise.setSize(s);damatrix.setSize(s);datuple.setSize(s);
				daderivative.setSize(s);dainverse.setSize(s);dainner.setSize(s);
					for(i=0,k=0;i<s;i++){
               //         coor[1] = i;
               //         for(coor[0] = 0 ; coor[0] < 1 ; coor[0]++) da_vectors(coor) = EM_obs[i].d.d[coor[0]];
						datuple[i] = EM_obs[i].d.d[0];
						inner_param[i].k = EM_obs[i].d.k;
				//		printf("%e\t%e\n",EM_obs[i].d.k[0],EM_obs[i].d.d[0]);
						for(j=0;j<i;j++,k++){
							double d = ExOp::pnorm(EM_obs[i].d.k[0] -EM_obs[j].d.k[0]);
							//damatrix.data[k] = exp(scope[1] - exp(scope[0]) * d*d);
							//daderivative.data[k] = damatrix.data[k] * (-exp(scope[0]) * d*d);
							damatrix_nonoise.data[k] =  exp(scope[1] - exp(scope[0]) * d);
							daderivative.data[k] = damatrix_nonoise.data[k] * (-d);
						}
						damatrix_nonoise.data[k] = exp(scope[1]);
						daderivative.data[k] = 0.0f;
						k++;
					}
					damatrix = damatrix_nonoise;
             //       for(i=0,j=0;j<s;j++) {damatrix.data[i] += exp(scope[2]+ scope[1]);i += j+2;}
                    for(i=0,j=0;j<s;j++) {damatrix.data[i] += exp(scope[2]);i += j+2;}
	if (showmat)	{Matrix<double>(damatrix).show(); printf("%e\t%e\n",scope[0],exp(scope[0]));}

		dainverse = damatrix.mkInverse();
	if (!ExOp::isValid(dainverse.data[0])) dainverse = damatrix.inverse_MK2();



	//Matrix<double>(damatrix_nonoise).show();


	//	(Matrix<double>(damatrix) * Matrix<double>(dainverse)).show();
                //    if (fabs(test.data[6] - 1.0f) > 0.001f) {damatrix.show(); printf("OMG!\n"); }
                    //test.show();

				//	datuple.show();
                 //   dainverse.show();
				//	(dainverse * datuple).show();
//	ExOp::show(datuple);
					dainner = dainverse.Xformed_outer_product(datuple);
                 //   printf("davec: ");
                 //   ExOp::show(dainverse * datuple );

                //    printf("damat: ");dainner.show();
					dainner -= dainverse;

				/*	dainner *= exp(-scope[1]);
					dainverse *= exp(-scope[1]);
					damatrix *= exp(scope[1]);
					damatrix_nonoise *= exp(scope[1]);
					daderivative *= exp(scope[1]);*/
				//	ExOp::show( dainverse.Xformed_inner_product(datuple) );
				//	ExOp::show(dainner.trace());
				//	ExOp::show(dainner.trace_of_product(damatrix));
					i_deriv[2] = 0.5f * exp(scope[2]) * dainner.trace(); // +scope[1] in exp
                    i_deriv[1] = 0.5f * dainner.trace_of_product(damatrix_nonoise) ; // + i_deriv[2]
					i_deriv[0] = 0.5f * exp(scope[0]) * dainner.trace_of_product(daderivative);
                //    printf("%e\t%e\t%e    .. %e   %e \n", deriv[0], deriv[1], deriv[2],dainverse.Xformed_inner_product(datuple),damatrix.log_determinant());

              //   ttmptmp = -0.5f * (dainverse.Xformed_inner_product(datuple) + damatrix.log_determinant());

               	  // ttmptmp = dainverse.Xformed_inner_product(datuple);
               	  ttmptmp = damatrix.Xformed_inner_product_of_inverse(datuple);


               	  if (ttmptmp < 0.0f) ExOp::toMin(ttmptmp);
               	  else ttmptmp = -0.5f * (damatrix.log_determinant() + ttmptmp);

	       //         printf("damat: %e\n", damatrix.log_determinant());
              // 	  ttmptmp = damatrix.log_determinant();
              //  if (ttmptmp > 500) printf("%e\t%e\n",dainverse.Xformed_inner_product(datuple), damatrix.log_determinant() ); //	(Matrix<double>(damatrix) * Matrix<double>(dainverse)).show();
	for(i=0;i<NBPARAM*2+1;i++){
		tmp_dev = (scope[i] -  param_prior_mean[i]) / param_prior_std[i];
		ttmptmp -= 0.5f * tmp_dev * tmp_dev; i_deriv[i] -= tmp_dev;
	}

	Tuple<double,0u> datuple2 = dainverse * datuple;

	for(i=0;i<s;i++) inner_param[i].d[0] = datuple2[i];

		return ttmptmp;
    }

LFHTEMP	double GaussianProcessProbSimple<NBPARAM,NBINPUT>::EMfinit(){


//	if (EMstep_count == 0){

//		llikelyhood = noise.EMfinit();
//		scope[1] = log(noise.getVar(0));
//		scope[2] = log(noise.getVar(0));
//	}
//		noise.EMfinit();
	// WHERE THE FUN HAPPENS !!!


        double eps = 0.000001f;

        llikelyhood = setDerivative_routine(deriv);

//	KeyElem< Tuple<double, NBPARAM> , Tuple<double, NBINPUT> > * inner_param;


  //     double fake[3];double tmp;
  //      for(unsigned int i =0;i<3;i++){ tmp = scope[i]; scope[i] += eps; printf("%i: Rel Err = %e  (%e)\n",i, 1.0f  - ((setDerivative_routine(fake)- llikelyhood) / (deriv[i] *eps) ), deriv[i] ); scope[i] = tmp;}





/*

					dainverse = damatrix.inverse();

				//	(Matrix<double>(damatrix) * Matrix<double>(dainverse) ).show();
					datuple.show();
					(dainverse * datuple).show();
					dainner = dainverse.Xformed_outer_product(datuple);
					dainner -= dainverse;

					dainner.show();	printf("%f\n", dainner.trace());

					deriv[0] = -0.5f * dainner.trace();
					deriv[1] = -0.5f * (dainner.trace_of_product(damatrix) - dainner.trace());
					deriv[2] = -0.5f * dainner.trace_of_product(daderivative);
					value = -0.5f * (dainverse.Xformed_inner_product(datuple) + damatrix.log_determinant());

					printf("deriv = [%e,%e,%e]\n", deriv[0]*0.000001f,deriv[1]*0.000001f,deriv[2]*0.000001f);
					printf("F[%e,%e,%e] = %e\n", guess[0],guess[1],guess[2], value);
*/

return llikelyhood;}

LFHTEMP	void GaussianProcessProbSimple<NBPARAM,NBINPUT>::EMnext(double alpha, GaussianProcessProbSimple<NBPARAM, NBINPUT> &fout) const{
     fout.alpha_search = alpha_search;
     memcpy(fout.scope, scope, sizeof(scope) );
     fout.alpha_search.updateAscent(llikelyhood, fout.scope, deriv);
}
LFHTEMP	void GaussianProcessProbSimple<NBPARAM,NBINPUT>::EMstep(double alpha){
  //   printf("step %e\n", alpha_search.updateAscent(llikelyhood, scope, deriv));
	printf("bef: %e\t%e\t%e\n", scope[0], scope[1], scope[2]);
    printf("%e \t", alpha_search.updateAscent(llikelyhood, scope, deriv));
	printf("aft: %e\t%e\t%e\n", scope[0], scope[1], scope[2]);
}

LFHTEMP	Tuple<double,NBINPUT> GaussianProcessProbSimple<NBPARAM,NBINPUT>::getSignalMean(const Tuple<double, NBPARAM>& whe) const{ Tuple<double,NBINPUT> fout;
    // needs to compute distance to all key points
    unsigned int i,j,k,l;
    Tuple<double, NBPARAM> dif;
    for(i=0 ; i < inner_param.getSize();i++){
        dif = whe - inner_param[i].k;
        // multiply my upper triangular and computes self-innerproduct:
        if (i == 0) fout = inner_param[0].d * exp(scope[1] - exp(scope[0]) * ExOp::pnorm(dif) );
        else fout +=  inner_param[i].d * exp(scope[1] - exp(scope[0]) * ExOp::pnorm(dif) );
    }
return(fout);}
LFHTEMP	void GaussianProcessProbSimple<NBPARAM,NBINPUT>::setGrid(const Vector< Tuple<double, NBPARAM> > &ingr){
	if (ingr.size() == 0) {static_warning_handdle << LFH_WARNING_UNEXPECTED_INPUT;return;}
	resample_to_fixgrid = true;
	nb_param = ingr.size();
	inner_param.setSize(nb_param);
	KeyElem< Tuple<double,NBPARAM>,  Tuple<double,NBINPUT> > a;
	ExOp::toZero(a);
	for(int i=0;i<ingr.size();i++) {a.k = ingr[i]; inner_param[i] = a;}
}

LFHTEMP	void GaussianProcessProbSimple<NBPARAM,NBINPUT>::show(FILE* f, int level)const{

}





/*

double LL_gradient(double x, double r, double p){}

#undef LFHTEMP
#define LFHTEMP template<class F, class A, class G>

LFHTEMP template<Tuple_flag TF> double InflatedNegativeBinomialFunction<F,A,G>::operator()(uint32_t value, const Tuple<Tuple<double, 0u, TF>, 2u> &hidden)const{
    double rval = mat_r(hidden);
    double pval = mat_p(hidden);


}*/

