/*
 * primitive.cpp
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

#include "primitive.h"

namespace LFHPrimitive{
	int verbose;
	int my_tmpvariable_for_error_msgs;
	myHashmap<uint32_t, void*> lfhstatic_ressources;
	// 8x12 font for 96 characters
	unsigned int def_font[] = {
	0x00000000,0x00000000,0x00000000,// space
	0x00180000,0x18181818,0x00181818,// !
	0x00000000,0x00000000,0x003C3C3C,
	0x12120000,0x2424247F,0x004848FE,
	0x3C180000,0x063C6066,0x00183C66,
	0xF3600000,0x30180C66,0x0006CF66,
	0x33FE0000,0x0C0EDB7B,0x001C0606,
	0x00000000,0x00000000,0x00181818,
	0x0C181830,0x0C0C0C0C,0x00301818,
	0x3018180C,0x30303030,0x000C1818,
	0x00000000,0x3CFF3C66,0x00000066,
	0x18180000,0x1818FF18,0x00000018,
	0x00180C00,0x00000000,0x00000000,
	0x00000000,0x0000FF00,0x00000000,
	0x00180000,0x00000000,0x00000000,
	0x03000000,0x30180C06,0x0000C060,

	0x663C0000,0x66666666,0x003C6666,
	0x187E0000,0x18181818,0x00181E18,
	0x067E0000,0x6030180C,0x003C6660,
	0x663C0000,0x60386060,0x003C6660,
	0x30780000,0x3C3C367E,0x00303838,
	0x663C0000,0x063E6060,0x007E0606,
	0x663C0000,0x063E6666,0x00380C06,
	0x0C0C0000,0x30301818,0x007E6660,
	0x663C0000,0x663C6666,0x003C6666,
	0x301C0000,0x667C6060,0x003C6666,
	0x00180000,0x18000000,0x00000000,
	0x00180C00,0x18000000,0x00000000,
	0x18300000,0x180C060C,0x00000030,
	0x00000000,0x7E007E00,0x00000000,
	0x180C0000,0x18306030,0x0000000C,
	0x00180000,0x60301818,0x003C6660,

	0x063C0000,0xA5A5BDDB,0x007CC6BB,
	0x66E70000,0x3C24667E,0x00101818,
	0xC67F0000,0xC67EC6C6,0x007FC6C6,
	0xC67C0000,0x03030303,0x007CC603,
	0x663F0000,0xC6C6C6C6,0x003F66C6,
	0xC6FF0000,0x363E3606,0x00FFC606,
	0x061F0000,0x363E3606,0x00FFC606,
	0xC67C0000,0x0303F3C3,0x007CC603,
	0x66E70000,0x667E6666,0x00E76666,
	0x187E0000,0x18181818,0x007E1818,
	0x331E0000,0x30303033,0x007C3030,
	0x66EF0000,0x361E1E36,0x00EF6636,
	0xCCFF0000,0x0C0C0C0C,0x003F0C0C,
	0x66E70000,0x7E5A5A42,0x00E7667E,
	0x666F0000,0x6E7E7676,0x00F7666E,
	0x663C0000,0xC3C3C3C3,0x003C66C3,

	0x061F0000,0xC67E0606,0x007FC6C6,
	0x663CFC00,0xC3C3C3C3,0x003C66C3,
	0x66EF0000,0xC67E3636,0x007FC6C6,
	0xC37E0000,0x037EC0C0,0x007EC303,
	0x183C0000,0x18181818,0x00FFDB18,
	0x663C0000,0x66666666,0x00FF6666,
	0x18180000,0x663C3C3C,0x00E76666,
	0x24240000,0x425A7E7E,0x00E76666,
	0x66E70000,0x3C183C3C,0x00E7663C,
	0x183C0000,0x3C3C1818,0x00E76666,
	0xC6FF0000,0x30180C0C,0x00FF6330,
	0x0C0C0C3C,0x0C0C0C0C,0x003C0C0C,
	0xC0000000,0x0C183060,0x00000306,
	0x3030303C,0x30303030,0x003C3030,
	0x00000000,0x00000000,0x183C6600,
	0x000000FF,0x00000000,0x00000000,

	0x00000000,0x00000000,0x000C1800,
	0x63FE0000,0x3E607E63,0x00000000, // a
	0xC67F0000,0x7EC6C6C6,0x00070606,
	0xC37E0000,0x7EC30303,0x00000000,
	0x63FE0000,0x7E636363,0x00706060,
	0xC37E0000,0x7EC3FF03,0x00000000, // e
	0x0C3E0000,0x3E0C0C0C,0x00380C0C,
	0x7E60603E,0xFE636363,0x00000000,
	0x66EF0000,0x3E6E6666,0x00070606,
	0x187E0000,0x1E181818,0x00180000,
	0x3030301E,0x3E303030,0x00300000,
	0x36E70000,0x76361E1E,0x00070606,
	0x187E0000,0x18181818,0x001C1818,
	0x4ADB0000,0x3F6A6A6A,0x00000000,
	0x66E70000,0x3F6E6666,0x00000000,
	0xC37E0000,0x7EC3C3C3,0x00000000,

	0xC67E060F,0x7FC6C6C6,0x00000000,
	0x637E60F0,0xFE636363,0x00000000,
	0x0C3F0000,0x7FDC0C0C,0x00000000,
	0xC37E0000,0x7EC31E70,0x00000000,
	0x6C380000,0x3E0C0C0C,0x00000C0C,
	0x76FC0000,0x77666666,0x00000000,
	0x3C180000,0xE766663C,0x00000000,
	0x3C240000,0xC3665A7E,0x00000000,
	0x66C30000,0xC3663C3C,0x00000000,
	0x3C18180E,0xC366663C,0x00000000,
	0x667F0000,0x7F33180C,0x00000000,
	0x18181870,0x18180E18,0x00701818,
	0x18181818,0x18181818,0x18181818,
	0x1818180E,0x18187018,0x000E1818,
	0x00000000,0x73000000,0x0000CEDB,
	0x00000000,0x00000000,0x00000000,

	0xC37E0000,0x7EC3FF03,0x60FC6000,
	0xC37E0000,0x7EC3FF03,0x00FC0000};


	myHashmap< unsigned int , pair< void*, unsigned int> > central_alias_bank;
	myHashmap< void* , unsigned int > central_alias_bank_back;

	void LFHBreakpoint(){

	}

	string to_string(int i){
		char buffer[16];
		sprintf(buffer, "%i",i);
		return string(buffer);
	}

	char* cloneString(const char* const what){
		unsigned int i = strlen(what) +1;
		char* _out = new char[i];
		memcpy(_out,what,sizeof(char)*i);
		return(_out);
	}
// exp(-x) + exp(-2x)+ exp(-2x))

// log(incbeta(a,b,1 / (1 + exp(-x))))

Vector<Tuple<double, 5u> > derivativeAnalyser(double (&fnc)(const double),double (&der)(const double), bool posonly, bool show){ Vector<Tuple<double, 5u> > fout;
    Tuple<double, 5u> fout_in;
    double values[5];
    int expon[4];
    double offset;
    fout.setSize(800);
    int32_t i = (posonly) ? 400 : 0;
    for(; i< 800; i++){

        fout_in[0] = (i < 400) ? -3.0 * pow(2.0, 200 - i) : 3.0 * pow(2.0, -600 + i);
        fout_in[1] = fnc(fout_in[0]); fout_in[2] = der(fout_in[0]);
        values[0] = frexp(fout_in[1], expon); values[1] = frexp(fout_in[2], expon+1);
        offset = pow(0.5, 30 - expon[0] + expon[1]) * values[0] / values[1]; // expected position for small but noticable change

        values[2] = fnc(fout_in[0]+offset); values[3] = fnc(fout_in[0]- offset);
        values[2] = frexp(values[2] - values[3], expon+2);
        fout_in[3] = (pow(0.5, 29 -expon[0] + expon[2]) * values[2] / values[0]) - 1.0;

        expon[3] = (i < 400) ? 170 - i :  -570 + i;
        offset = pow(2.0, expon[3]);
        values[2] = fnc(fout_in[0]+offset); values[3] = fnc(fout_in[0]- offset); // expon[3] * values[1]
        values[2] = frexp(values[2] - values[3], expon+2);
        fout_in[4] = (pow(0.5, -1-expon[1] + expon[2] - expon[3]) * values[2] / values[1]) - 1.0;



        // F[x] = F[0] + F'[0] x
// K = (F[x] - F[0]) / F[0]  = F'[0] x / F[0]    =>   x  = K * F[0] / F'[0]
        fout[i] = fout_in;

    }
    if (show){
        i = (posonly) ? 400 : 0;
        for(i=0;i<fout.getSize();i++) printf("F/F'[%e]->%e,%e\tRelErr-1:%e\t%e\n", fout[i][0], fout[i][1], fout[i][2], fout[i][3], fout[i][4]);
    }
return fout;}

double deflatedNBterm_approx(const double r, const double m){ // r - log(exp(m)+1) -log(exp(log(1+exp(-m))*exp(r)) -1.0f);
   // if (x >= -3.0) return x - log(x + log(1.0+exp(-x)));
    return  r- log(exp(m)+1) -log(exp(exp(quasilimit_func(m)-m+r)) -1.0f);
}
// exp(x + log(log(1.0+exp(-x)))) = exp(m) * log(1.0+exp(-m))

// r - log(exp(m)+1) -log(exp( log(1+exp(-m)) *exp(r)) -1.0f);
// or if log(1+exp(-m)) *exp(r) <<< epsilon
// tmp = log(1+exp(-m)) *exp(r)
// r - log(exp(m)+1) -log( tmp + tmp^2/2);
// - log(exp(m)+1) -log( (log(1+exp(m))-m) * (1.0 + 0.5 * log(1+exp(-m)) *exp(r) );


double deflatedNBterm(const double r, const double m){ // r - log(exp(m)+1) -log(exp( log(1+exp(-m)) *exp(r) ) -1.0f);
    /*if (m > 18.375f) return r-m -log(exp( exp(r-m) ) -1.0f);
    if (m - r > 36.75f){
        if (m > 18.375f) return 0;
        else return -log(exp(m)+1) -log(log(1+exp(-m)));
    }
else*/
    //if (m - r > 18.375f) return r - log(exp(m)+1) -log(exp(log(1+exp(-m))*exp(r)) -1.0f);
    /*if (r < -18.375f){
        if (m <= 3.0) return -log(exp(m)+1) -log(log(1+exp(m))-m);
        else {
            return  -log(exp(-m)+1) -quasilimit_func(-m);
        }
    }


    if (m > 18.375f){
        return r-m -log(exp( exp(r-m) ) -1.0f);
    }
    double tmp = quasilimit_func(-m) + r-m;
    if (tmp >= 18.375f){
        return r - log(exp(m)+1) - exp(tmp) - exp(-exp(tmp));

    }*/
    //return r - log(exp(m)+1) -log(exp(exp(tmp)) -1.0f);
    //double tmp = (x <= 3.0)
    //else
   // double s;
    if (m <= 3){
        if (r >= -18.375f) return r - log(exp(m)+1) -log_e_x_m1((log(exp(m)+1.0)-m) *exp(r));
      //  s = (log(exp(m)+1.0)-m) *exp(r);
      //  return (m-log(exp(m)+1.0)) * (1.0 + exp(r) * (0.5 + s * (1.0/24.0 - s * s *(1.0/2880.0 + s * (1.0/840 - s * (217.0/508032 -s * (635229.0/5080320 - s* (5039857.0/29030400))))))) - m - log(log(exp(m)+1.0)-m));
        return  - log(log(exp(m)+1.0)-m) -m  + (m-log(exp(m)+1.0)) *(1.0 + 0.5 * exp(r));
    }
    double t = exp(-m);
    t *= (0.5L - t * (5.0L/24.0L - t*(0.125L - t* (251.0L/2880.0L - t*(19.0L/288.0L - t * (19087.0L/362880.0L - t * (751.0L/17280.0L - t * (1070017.0L/29030400.0L - t * (2857.0L/89600.0L)))))))));
    if (r - m - t >= -1.0) return r - m - log(exp(-m)+1) -log_e_x_m1(exp(r-m-t));
   /*(r - m < -18.375f) ? -log(exp(-m)+1.0) -log_1minus_e_x(exp(-t)) :*/
    //return -log(exp(-m)+1.0)  log(exp(exp(r-m-t))-1.0)
    //t = exp(r-m-t) -r+m;
    //return  -log(exp( exp(n) * exp(-t) + p) - exp(p));

    /*1  + exp(n) * exp(-t) + p + (exp(n) * exp(-t) + p)^2
    1  + p + p^2 +*/

    //return 1.0f;
    //log(x+x^2/2 ...)

    //x = exp(x)
  //  return ((m > 18.375) ? -exp(-m): -log(exp(-m)+1)) + r-m-log_1minus_e_x(exp(r-m-t));
    double s = exp(r-m-t);
    return ((m > 18.375) ? -exp(-m): -log(exp(-m)+1)) + t - s * (0.5 + s * (1.0/24.0 - s * s *(1.0/2880.0 + s * (1.0/840 - s * (217.0/508032 -s * (635229.0/5080320 - s* (5039857.0/29030400)))))));



    //log(exp(hehe+m-r) - exp(m-r))
    // return -quasilimit_func(-m) + r-m;
    //  -log(exp( log(1+exp(-m)) *exp(r) ) -1.0f);
} // log(exp(x) -1) x + log(1.0 - exp(-x)) +- x - exp(-x)
// exp(-quasilimit_func(x)) = exp(-x)*log(1.0+exp(x))
double quasilimit_func(const double x){ // x - log(log(1.0+exp(x)))
    if (x >= -3.0) return x - log(x + log(1.0+exp(-x)));
    double t = exp(x);
    return t * (0.5 - t * (5.0/24.0 - t*(0.125 - t* (251.0/2880.0 - t*(19.0/288.0 - t * (19087.0/362880.0 - t * (751.0/17280.0 - t * (1070017.0/29030400.0 - t * (2857.0/89600.0)))))))));
}
double quasilimit_func_dx(const double x){ // x + log(log(1.0+exp(-x)))

    if (x >= -3.0) return 1.0 - 1.0 / ((x + log(1.0+exp(-x)))*(1.0+exp(-x)));
    double t = exp(x);
    return t * (0.5 - t * (5.0/12.0 - t*(0.375 - t* (251.0/720.0 - t*(95.0/288.0 - t * (6*19087.0/362880.0 - t * (7*751.0/17280.0 - t * (8*1070017.0/29030400.0 - t * (9*2857.0/89600.0)))))))));
}
double deflatedNBterm_dr(const double r, const double m){
    /*if (m <= 3){
        //if (r < -18.375) return -exp(-r) / ((log(exp(m)+1.0)-m) + exp(-r));
        //else
        return -one_over_e_x_1minus((log(exp(m)+1.0)-m) *exp(r));
    }*/
    double s;
    if (m <= 3){
        if (r >= -18.375) return 1.0 - (log(exp(m)+1.0)-m) * exp(r) / (1.0 - exp(-(log(exp(m)+1.0)-m) *exp(r)));
        s = (log(exp(m)+1.0)-m) *exp(r);
        return -s * (0.5 + s * (2.0/24.0 - s * s *(4.0/2880.0 + s * (5.0/840 - s * (1302.0/508032 -s * (635229.0*7/5080320 - s* (8*5039857.0/29030400)) )))));
        //return -log(exp(m)+1.0) -log(log(exp(m)+1.0)-m);
    }

    double t = exp(-m);
    t *= (0.5L - t * (5.0L/24.0L - t*(0.125L - t* (251.0L/2880.0L - t*(19.0L/288.0L - t * (19087.0L/362880.0L - t * (751.0L/17280.0L - t * (1070017.0L/29030400.0L - t * (2857.0L/89600.0L)))))))));
    if (r - m - t >= -1.0) return 1.0 - exp(r-m-t) / (1.0 - exp(-exp(r-m-t)));
    s = exp(r-m-t);
    return -s * (0.5 + s * (2.0/24.0 - s * s *(4.0/2880.0 + s * (5.0/840 - s * (1302.0/508032 -s * (635229.0*7/5080320 - s* (8*5039857.0/29030400)) )))));
/*
//-1.0 / (exp(tmp) - 1.0);

    // (-exp(-log(1+exp(-m))*exp(r)) -(log(1+exp(-m))*exp(r)+1.0f))/(1.0f -exp(log(1+exp(-m))*exp(r)));
    if (r < -18.375f){
        return (-exp(-log(1+exp(-m))*exp(r)) -log(1+exp(-m))*exp(r)+1.0f) / (1.0f -exp(log(1+exp(-m))*exp(r)));
    }

    double tmp = log(1+exp(-m))*exp(r);
    if (r > 18.375) return exp(r-m);

    if (tmp >= 18.375f){
        return 1.0 - tmp / (1.0f - exp(-tmp));
    }
    return 1.0 - exp(tmp) * tmp / (exp(tmp) -1.0);*/
}
double deflatedNBterm_dm(const double r, const double m){ // - 1/(1 + e^-m) -log(exp(log(1+e^-m)e^r) -1)

    if (m <= 3){
        if (r < -18.375) return (1.0 - exp(m)*(log(exp(m)+1)-m)) / ((exp(m)+1)*(log(exp(m)+1)-m));
        //return -1.0 / (exp(-m)+1) - ( log(exp(-m)+1.0) *exp(r) ) / ( log(exp(-m)+1.0) *exp(r) );
        //return -exp(m) / (exp(m)+1) + (log(exp(-m)+1.0)) *exp(r) / ((exp(m)+1)*(1.0 - exp(-(log(exp(m)+1.0)-m) *exp(r))));
        //return -exp(m) / (exp(m)+1) + exp((log(exp(m)+1.0)-m) *exp(r))*(log(exp(m)+1.0)-m) *exp(r) / ((1+exp(m))*(exp((log(exp(m)+1.0)-m) *exp(r)) -1.0));
        return -exp(m) / (exp(m)+1) + exp(r) / ((exp(m)+1)*(1.0 - exp(-(log(exp(m)+1.0)-m) *exp(r)) ));
        //return (1.0 - exp(m)*log(1+exp(-m))) / ((exp(m)+1)*(log(1+exp(m))-m));
    } // return -log(exp(((log(exp(m)+1.0)-m) *exp(r))) - 1.0);

    double t = exp(-m);
    double dt = t * (0.5 - t * (5.0/12.0 - t*(0.375 - t* (251.0/720.0 - t*(95.0/288.0 - t * (6.0*19087.0/362880.0 - t * (7.0*751.0/17280.0 - t * (8.0*1070017.0/29030400.0 - t * (9.0*2857.0/89600.0)))))))));
    t *= (0.5L - t * (5.0L/12.0L - t*(0.125L - t* (251.0L/2880.0L - t*(19.0L/288.0L - t * (19087.0L/362880.0L - t * (751.0L/17280.0L - t * (1070017.0L/29030400.0L - t * (2857.0L/89600.0L)))))))));

    if (r - m - t >= 0.0l) return - 1.0 - 1.0/(exp(-m)+1) - exp(r-m-t) * (dt-1)/ (1.0 - exp(-exp(r-m-t)));

    double s = exp(r-m-t);
    return 1.0 / (1.0 + exp(m)) - dt + s * (0.5 + s * (2.0/24.0 - s * s *(4.0/2880.0 + s * (5.0/840 - s * (1302.0/508032 -s * (635229.0*7/5080320 - s* (8*5039857.0/29030400)) )))));

}
double deflatedNBterm_dr2(const double r, const double m){ // - 1/(1 + e^-m) -log(exp(log(1+e^-m)e^r) -1)
    double s;
    if (m <= 3){
        if (r >= -12.375){
            s = exp(r) /(1.0 - exp(-(log(exp(m)+1.0)-m) *exp(r)));
            return - (log(exp(m)+1.0)-m) *(s -  s * s * exp(-(log(exp(m)+1.0)-m) *exp(r)) * (log(exp(m)+1.0)-m) );

        }
//            return              (log(exp(m)+1.0)-m) * exp(r) * (1.0 +
//+ exp(-(log(exp(m)+1.0)-m) *exp(r)) * ((log(exp(m)+1.0)-m) * exp(r) - 1.0) ) / ((-1.0 + exp(-(log(exp(m)+1.0)-m) *exp(r))) * (1.0 - exp(-(log(exp(m)+1.0)-m) *exp(r)))); // TOCHECK
        s = -(log(exp(m)+1.0)-m) *exp(r);
        return s * (0.5 + s * (4.0/24.0 - s * s *(16.0/2880.0 + s * (25.0/840 - s * (6.0*1302.0/508032 -s * (635229.0*49.0/5080320 - s* (5039857.0*64.0/29030400)) )))));
        //return -log(exp(m)+1.0) -log(log(exp(m)+1.0)-m);
    }

    double t = exp(-m);
    t *= (0.5L - t * (5.0L/24.0L - t*(0.125L - t* (251.0L/2880.0L - t*(19.0L/288.0L - t * (19087.0L/362880.0L - t * (751.0L/17280.0L - t * (1070017.0L/29030400.0L - t * (2857.0L/89600.0L)))))))));
    if (7.0f +r >= m+t) return (-exp(r-m-t) / (1.0 - exp(-exp(r-m-t))) + exp(-exp(r-m-t)) * exp(2.0 * (r-m-t)) / ((1.0 - exp(-exp(r-m-t)))*(1.0 - exp(-exp(r-m-t))))); // TOCHECK
    s = exp(r-m-t);
    return -s * (0.5 + s * (4.0/24.0 - s * s *(16.0/2880.0 + s * (25.0/840 - s * (6.0*1302.0/508032 -s * (635229.0*49.0/5080320 - s* (5039857.0*64.0/29030400)) )))));
}
double deflatedNBterm_dm2(const double r, const double m){ // - 1/(1 + e^-m) -log(exp(log(1+e^-m)e^r) -1)
    double s;
    if (m <= 3){
        if (r < -12.0f){ //return (1.0 - exp(m)*(log(exp(m)+1)-m)) / ((exp(m)+1)*(log(exp(m)+1)-m));
            s = (log(exp(m)+1)-m);
            return (1.0 - exp(m) * s *(s + 1.0)) / (s * s * (exp(m) +1) * (exp(m) +1));
        }    //             return  exp(m)*(( ((exp(m)+1)*(log(exp(m)+1)-m))*(exp(m) / (exp(m)+1) - (log(exp(m)+1) - m) )) - exp(m)*(1.0 - exp(m)*(log(exp(m)+1))-m) * ((log(exp(m)+1)-m) + (exp(m)+1)/(exp(m)+1) ) )/ (((exp(m)+1)*(log(exp(m)+1)-m))*((exp(m)+1)*(log(exp(m)+1)-m)));
        //return -1.0 / (exp(-m)+1) - ( log(exp(-m)+1.0) *exp(r) ) / ( log(exp(-m)+1.0) *exp(r) );
        //return -exp(m) / (exp(m)+1) + (log(exp(-m)+1.0)) *exp(r) / ((exp(m)+1)*(1.0 - exp(-(log(exp(m)+1.0)-m) *exp(r))));
        //return -exp(m) / (exp(m)+1) + exp((log(exp(m)+1.0)-m) *exp(r))*(log(exp(m)+1.0)-m) *exp(r) / ((1+exp(m))*(exp((log(exp(m)+1.0)-m) *exp(r)) -1.0));
        s =exp(-(log(exp(m)+1.0)-m) *exp(r));
        return exp(m) * (-exp(r)-1.0 - s * (s - exp(2.0*r-m) - 2.0 - exp(r) )) / ((s-1)*(s-1)*(exp(m)+1)*(exp(m)+1));
        //return exp(m) * (-1.0 - s*(s*(exp(m)+1) - exp(2.0*r-m) - 2.0 - exp(r))) / ((s-1)*(s-1)*(exp(m)+1)*(exp(m)+1));
        //return ((exp(r*2+m)* exp(-(log(exp(m)+1.0)-m) *exp(r))/ (1.0 + exp(m))) / ((1.0 - exp(-(log(exp(m)+1.0)-m) *exp(r)))*(1.0 - exp(-(log(exp(m)+1.0)-m) *exp(r)))) -exp(m))* (exp(m)+1) -(s -exp(m)) / ((exp(m)+1)*(exp(m)+1));
        //return (1.0 - exp(m)*log(1+exp(-m))) / ((exp(m)+1)*(log(1+exp(m))-m));
    } // return -log(exp(((log(exp(m)+1.0)-m) *exp(r))) - 1.0);

    double t = exp(-m);
    double dt = t * (0.5 - t * (5.0/12.0 - t*(0.375 - t* (251.0/720.0 - t*(95.0/288.0 - t * (6.0*19087.0/362880.0 - t * (7.0*751.0/17280.0 - t * (8.0*1070017.0/29030400.0 - t * (9.0*2857.0/89600.0)))))))));//TODO
    double ddt = t * (0.5 - t * (10.0/12.0 - t*(0.375*3.0 - t* (251.0*4.0/720.0 - t*(95.0*5.0/288.0 - t * (36.0*19087.0/362880.0 - t * (49.0*751.0/17280.0 - t * (64.0*1070017.0/29030400.0 - t * (81.0*2857.0/89600.0)))))))));//TODO
    t *= (0.5L - t * (5.0L/12.0L - t*(0.125L - t* (251.0L/2880.0L - t*(19.0L/288.0L - t * (19087.0L/362880.0L - t * (751.0L/17280.0L - t * (1070017.0L/29030400.0L - t * (2857.0L/89600.0L)))))))));//TODO

    if (r - m - t >= 0.0l) return  1.0 + 1.0/(exp(-m)+1) + exp(r-m-t) * (ddt + (dt-1))/ (1.0 - exp(-exp(r-m-t))) - exp(-exp(r-m-t)) * exp(2*(r-m-t)) * (dt-1)/ ((1.0 - exp(-exp(r-m-t)))*(1.0 - exp(-exp(r-m-t)))); //TODO

     s = exp(r-m-t);
    return -1.0 / (1.0 + exp(m)) + dt - s * (0.5 + s * (4.0/24.0 - s * s *(16.0/2880.0 + s * (25.0/840 - s * (6.0*1302.0/508032 -s * (635229.0*49.0/5080320 - s* (64.0*5039857.0/29030400)) )))));
}

void deflatedNBterm_wr_f_dr_dm(double* fout3, const uint32_t x, const double r, const double m){ // !x ignored!

    double s,t;
    if (m <= 3.0){
        if (r >= -12.0f){
            fout3[0] = r - log(exp(m)+1) -log_e_x_m1((log(exp(m)+1.0)-m) *exp(r));
            s = exp(r) /(1.0 - exp(-(log(exp(m)+1.0)-m) *exp(r)));
            fout3[1] = 1.0 - (log(exp(m)+1.0)-m) * s;
            fout3[2] = (s-exp(m)) / (exp(m)+1);
        }else{
            s = (log(exp(m)+1.0)-m) *exp(r);
            fout3[0] = -log(log(exp(m)+1.0)-m) -log(exp(m)+1.0) - 0.5 * s;
            fout3[1] = -s * (0.5 + s * (2.0/24.0 - s * s *(4.0/2880.0 + s * (5.0/840 - s * (1302.0/508032 -s * (635229.0*7.0/5080320 - s* (8.0*5039857.0/29030400)) )))));
            fout3[2] = (1.0 - exp(m)*(log(exp(m)+1.0)-m)) / ((exp(m)+1)*(log(exp(m)+1.0)-m));
        }
        if (x == 1u) return;
        fout3[0] -= log(exp(m)+1)* (x - 1u);
        fout3[2] -= (exp(m)* (x - 1u))/(exp(m)+1);
    }else{
        t = exp(-m);
        fout3[2] = t * (0.5 - t * (5.0/12.0 - t*(0.375 - t* (251.0/720.0 - t*(95.0/288.0 - t * (19087.0/60480.0 - t * (7.0*751.0/17280.0 - t * (8.0*1070017.0/29030400.0 - t * (9.0*2857.0/89600.0)))))))));
        t *= (0.5L - t * (5.0L/24.0L - t*(0.125L - t* (251.0L/2880.0L - t*(19.0L/288.0L - t * (19087.0L/362880.0L - t * (751.0L/17280.0L - t * (1070017.0L/29030400.0L - t * (2857.0L/89600.0L)))))))));
        if (7.0f +r >= m+t){ // if (r - m - t >= 0.0l)
            t = exp(r-m-t);
            fout3[0] = r - m - log(exp(-m)+1) -log_e_x_m1(t);
            fout3[1] = -t / (1.0 - exp(-t));
            fout3[2] = -1.0 / (exp(-m)+1) + (fout3[2]-1.0) * fout3[1];
            fout3[1] += 1.0;
        }else{
            s = exp(r-m-t);
            fout3[0] = ((m > 18.375) ? -exp(-m): -log(exp(-m)+1)) + t - s * (0.5 + s * (1.0/24.0 - s * s *(1.0/2880.0 + s * (1.0/840 - s * (217.0/508032 -s * (635229.0/5080320 - s* (5039857.0/29030400)))))));
            fout3[1] = -s * (0.5 + s * (2.0/24.0 - s * s *(4.0/2880.0 + s * (5.0/840 - s * (1302.0/508032 -s * (635229.0*7/5080320 - s* (8*5039857.0/29030400)) )))));
            fout3[2] = 1.0 / (1.0 + exp(m)) - fout3[2] - fout3[1];
        }
        if (x == 1u) return;
        fout3[0] -= (m + log(exp(-m)+1)) * (x - 1u);
        fout3[2] -= ((double)(x - 1u)) /(exp(-m)+1);
    }


    if (x < 16u){
        if (r >= 18.375f){
            s = exp(-r);
            t = 2.0 * s + 1.0f;
            for(uint32_t i=1;i<x;i++, t+= s) {fout3[0] += log(t); fout3[1] += 1.0 / t;}
            fout3[0] += r * (x-1u);
        }else{
            t = 2.0 + exp(r);
            s = exp(r);
            for(uint32_t i=1;i<x;i++, t+= 1.0) {fout3[0] += log(t); fout3[1] += s / t;}
        }
        return;
    }
    t = exp(r);
    fout3[0] += lngamma(t + x) - lngamma(t + 1.0);
    fout3[1] += (d_lngamma_dx(t + x) - d_lngamma_dx(t + 1.0)) * t;
}



void deflatedNBterm_wr_f_d2(double* fout6, const uint32_t x, const double r, const double m){ // !x ignored!
    double s,t;
    double mterm = log_e_x_p1(m);
    double mterminv = log_e_x_p1(-m);

    if (m <= 6.0){
        if (r >= -12.0f){
            fout6[0] = r - mterm -log_e_x_m1(mterminv *exp(r));
            s = exp(r) /(1.0 - exp(-mterminv *exp(r)));
            fout6[1] = 1.0 - mterminv * s;
            fout6[2] = (s-exp(m)) / (exp(m)+1);
            fout6[3] = s * mterminv *(s * exp(-mterminv *exp(r)) * mterminv - 1.0);
            t = exp(-mterminv *exp(r));
            fout6[4] = exp(m) * (-exp(r)-1.0 - t * (t - exp(2.0*r-m) - 2.0 - exp(r) )) / ((t-1)*(t-1)*(exp(m)+1)*(exp(m)+1));
            fout6[5] = (s / (exp(m)+1.0)) * (1.0 - mterminv * exp(r) / ( (-1.0 + exp(mterminv *exp(r)))));
            //fout6[5] = (s / (exp(m)+1)) * s * (1.0 + mterminv / (-1.0 + exp((log(exp(m)+1.0)-m) *exp(r))));
        }else{
            t = mterminv;
            s = t *exp(r);
            fout6[0] = -log(mterminv) - mterm - 0.5 * s;
            fout6[1] = -s * (0.5 + s * (2.0/24.0 - s * s *(4.0/2880.0 + s * (5.0/840 - s * (1302.0/508032 -s * (635229.0*7.0/5080320 - s* (8.0*5039857.0/29030400)) )))));
            fout6[2] = (1.0 - exp(m) * t) / ((exp(m)+1)*t);
            fout6[3] = -s * (0.5 - s * (4.0/24.0 - s * s *(16.0/2880.0 - s * (25.0/840 + s * (6.0*1302.0/508032 +s * (635229.0*49.0/5080320 + s* (5039857.0*64.0/29030400)) )))));
            fout6[4] = (1.0 - exp(m) * t *(t + 1.0)) / (t * t * (exp(m) +1) * (exp(m) +1));
            fout6[5] = -fout6[3] / (t * (exp(m)+1));
        }
        if (x == 1u) return;
        fout6[0] -= mterm * (x - 1u);
        s = (exp(m)+1);
        fout6[2] -= (exp(m)* (x - 1u))/s;
        fout6[4] -= (exp(m)* (x - 1u))/(s*s);
    }else{
/*
    double ddt = //TODO
    t *= (0.5L - t * (5.0L/12.0L - t*(0.125L - t* (251.0L/2880.0L - t*(19.0L/288.0L - t * (19087.0L/362880.0L - t * (751.0L/17280.0L - t * (1070017.0L/29030400.0L - t * (2857.0L/89600.0L)))))))));//TODO

    if (r - m - t >= 0.0l) return  1.0 + 1.0/(exp(-m)+1) + exp(r-m-t) * (ddt + (dt-1))/ (1.0 - exp(-exp(r-m-t))) - exp(-exp(r-m-t)) * exp(2*(r-m-t)) * (dt-1)/ ((1.0 - exp(-exp(r-m-t)))*(1.0 - exp(-exp(r-m-t)))); //TODO

     s = exp(r-m-t);
    return -1.0 / (1.0 + exp(m)) + dt - s * (0.5 + s * (4.0/24.0 - s * s *(16.0/2880.0 + s * (25.0/840 - s * (6.0*1302.0/508032 -s * (635229.0*49.0/5080320 - s* (64.0*5039857.0/29030400)) )))));
*/
/*
   double t = exp(-m);
    double dt = t * (0.5 - t * (5.0/12.0 - t*(0.375 - t* (251.0/720.0 - t*(95.0/288.0 - t * (6.0*19087.0/362880.0 - t * (7.0*751.0/17280.0 - t * (8.0*1070017.0/29030400.0 - t * (9.0*2857.0/89600.0)))))))));//TODO
    double ddt = t * (0.5 - t * (10.0/12.0 - t*(0.375*3.0 - t* (251.0*4.0/720.0 - t*(95.0*5.0/288.0 - t * (36.0*19087.0/362880.0 - t * (49.0*751.0/17280.0 - t * (64.0*1070017.0/29030400.0 - t * (81.0*2857.0/89600.0)))))))));//TODO
    t *= (0.5L - t * (5.0L/12.0L - t*(0.125L - t* (251.0L/2880.0L - t*(19.0L/288.0L - t * (19087.0L/362880.0L - t * (751.0L/17280.0L - t * (1070017.0L/29030400.0L - t * (2857.0L/89600.0L)))))))));//TODO

    if (r - m - t >= 0.0l) return  1.0 + 1.0/(exp(-m)+1) + exp(r-m-t) * (ddt + (dt-1))/ (1.0 - exp(-exp(r-m-t))) - exp(-exp(r-m-t)) * exp(2*(r-m-t)) * (dt-1)/ ((1.0 - exp(-exp(r-m-t)))*(1.0 - exp(-exp(r-m-t)))); //TODO

     s = exp(r-m-t);
    return -1.0 / (1.0 + exp(m)) + dt - s * (0.5 + s * (4.0/24.0 - s * s *(16.0/2880.0 + s * (25.0/840 - s * (6.0*1302.0/508032 -s * (635229.0*49.0/5080320 - s* (64.0*5039857.0/29030400)) )))));
*/

        t = exp(-m);
        s = t * (0.5L - t * (5.0L/24.0L - t*(0.125L - t* (251.0L/2880.0L - t*(19.0L/288.0L - t * (19087.0L/362880.0L - t * (751.0L/17280.0L - t * (1070017.0L/29030400.0L - t * (2857.0L/89600.0L)))))))));
        fout6[2] = t * (0.5 - t * (5.0/12.0 - t*(0.375 - t* (251.0/720.0 - t*(95.0/288.0 - t * (19087.0/60480.0 - t * (7.0*751.0/17280.0 - t * (8.0*1070017.0/29030400.0 - t * (9.0*2857.0/89600.0)))))))));
        if ((r >= m+s)){ // if (r - m - t >= 0.0l)
            fout6[4] = t * (0.5 - t * (10.0/12.0 - t*(0.375*3.0 - t* (251.0*4.0/720.0 - t*(95.0*5.0/288.0 - t * (36.0*19087.0/362880.0 - t * (49.0*751.0/17280.0 - t * (64.0*1070017.0/29030400.0 - t * (81.0*2857.0/89600.0)))))))));
            t = exp(r-m-s);
            fout6[0] = r - mterm -log_e_x_m1(t);
            fout6[1] = -t / (1.0 - exp(-t));
            //fout6[4] = 1.0 + 1.0/(exp(-m)+1) + t * (fout6[4] + (fout6[2]-1))/ (1.0 - exp(-t)) - exp(-t) * t * t * (fout6[2]-1)/ ((1.0 - exp(-t))*(1.0 - exp(-t)));
            //-1.0 / (1.0 + exp(m)) + fout6[2] - t * (0.5 + t * (4.0/24.0 - t * t *(16.0/2880.0 + t * (25.0/840 - t * (6.0*1302.0/508032 -t * (635229.0*49.0/5080320 - t* (64.0*5039857.0/29030400)) )))));
             fout6[3] = fout6[1] + exp(-t) * t * t / ((1.0 - exp(-t))*(1.0 - exp(-t)));
            fout6[5] = (fout6[1] + exp(-t) * t * t / ((1.0 - exp(-t))*(1.0 - exp(-t)))) * (-1.0 - fout6[2]);
            fout6[4] = -exp(-m) / ((exp(-m)+1)*(exp(-m)+1)) - fout6[4] * fout6[1] + (fout6[2]-1.0) * fout6[5];
            fout6[2] = -1.0 / (exp(-m)+1) + (fout6[2]-1.0) * fout6[1];
            fout6[1] += 1.0;
               //        -exp(r-m-t) / (1.0 - exp(-exp(r-m-t))) + exp(-exp(r-m-t)) * exp(2.0 * (r-m-t)) / ((1.0 - exp(-exp(r-m-t)))*(1.0 - exp(-exp(r-m-t)))); // TOCHECK
        }else{
            fout6[0] = s;
            s = exp(r-m-s);
            fout6[0] += -mterminv - s * (0.5 + s * (1.0/24.0 - s * s *(1.0/2880.0 + s * (1.0/840 - s * (217.0/508032 -s * (635229.0/5080320 - s* (5039857.0/29030400)))))));
            fout6[1] = -s * (0.5 + s * (2.0/24.0 - s * s *(4.0/2880.0 + s * (5.0/840 - s * (1302.0/508032 -s * (635229.0*7/5080320 - s* (8*5039857.0/29030400)) )))));
            fout6[4] = -1.0 / (1.0 + exp(m)) + fout6[2] - s * (0.5 + s * (4.0/24.0 - s * s *(16.0/2880.0 + s * (25.0/840 - s * (6.0*1302.0/508032 -s * (635229.0*49.0/5080320 - s* (64.0*5039857.0/29030400)) )))));
            //1.0 + 1.0/(exp(-m)+1) + s * (fout6[4] + (fout6[2]-1))/ (1.0 - exp(-s)) - exp(-s) * s * s * (fout6[2]-1)/ ((1.0 - exp(-s))*(1.0 - exp(-s)));
            fout6[3] = -s * (0.5 + s * (4.0/24.0 - s * s *(16.0/2880.0 + s * (25.0/840 - s * (6.0*1302.0/508032 -s * (635229.0*49.0/5080320 - s* (5039857.0*64.0/29030400)) )))));
            fout6[5] = -fout6[3] * (1.0 + fout6[2]);
            fout6[2] = 1.0 / (1.0 + exp(m)) - fout6[2] - fout6[1];
        }
        if (x == 1u) return;
        fout6[0] -= mterm * (x - 1u);
        s = (exp(-m)+1);
        fout6[2] -= ((double)(x - 1u)) / s;
        fout6[4] -= ((double)(x - 1u)) * exp(-m) /(s*s);
    }
    if (x < 16u){

        if (r >= 18.375f){
            s = exp(-r);
            t = s + 1.0;
            for(uint32_t i=1;i<x;i++, t+= s) {fout6[0] += log(t); fout6[1] += 1.0 / t; fout6[3] += s * i / (t * t);  }
            fout6[0] += r * (x-1u);
        }else{
            s = exp(r);
            t = 1.0 + s;// F(r) = sum  exp(r) *k / (k + exp(r))*(k + exp(r))
            for(uint32_t i=1;i<x;i++, t+= 1.0) {fout6[0] += log(t); fout6[1] += s / t; fout6[3] += s * i / (t * t);}
        }
        return;
    }
    if (r >= 12.0f){ // assuming exp(r) >> x then:
        // assuming cur >> k >> 1
        // lngamma(cur+k) - lngamma(cur+1) =+-=  0.5* (1.0 -k) * ( (1.0 -k) /cur - 1.0/(cur+k) + 1.0/(6*(cur+1)*(cur+k))) + (k-1)*log(cur)  6t + 6x*t + 12xt^2
        t = exp(-r);
        fout6[0] += (-1.0+x) * r + 0.5* (1.0 -x) * ( (1.0 -x) * t - t/(1.0+x*t) + 0.0*t*t /(6*(1 + t)*(1+t*x)));
        fout6[1] += (-1.0+x) - 0.5 * (1.0 -x) * ( (1.0 -x) * t - t/(1.0+x*t) + x*t/((1.0+x*t)*(1.0+x*t)) + 0.0*2.0*t*t /(6*(1 + t)*(1+t*x)) +  (1.0 + x + 2.0*x*t) *t* t*t /(6*(1 + t)*(1+t*x)*(1 + t)*(1+t*x)) );
        fout6[3] += 0.5* (1.0 -x) * ( (1.0 -x) * t - t/(1.0+x*t) + 2.0*x*t/((1.0+x*t)*(1.0+x*t)) - 2.0*x*x*t*t/((1.0+x*t)*(1.0+x*t)*(1.0+x*t)) + 0.0*4*t*t /(6*(1 + t)*(1+t*x)) ); // not exact... oh well!
    }else{
        t = exp(r);
        fout6[0] += lngamma(t + x) - lngamma(t + 1.0);
        s = (d_lngamma_dx(t + x) - d_lngamma_dx(t + 1.0));
        fout6[1] += s * t;
        fout6[3] += (s + (d2_lngamma_dx2(t + x) - d2_lngamma_dx2(t + 1.0)) * t) * t;
    }
}


double deflatedNB_pvalue_silly(const uint32_t x, const double r, const double m){
    double s,t;
    double lead, fout, tmptmp;
    double mterm = log_e_x_p1(m);
    double mterminv = log_e_x_p1(-m);
    if (m <= 3.0){
        if (r >= -12.0f){
            lead = r - mterm -log_e_x_m1(mterminv *exp(r));
        }else{
            t = mterminv;
            s = t *exp(r);
            lead = -log(mterminv) -log(exp(m)+1.0) - 0.5 * s;
        }
    }else{
        t = exp(-m);
        s = t * (0.5L - t * (5.0L/24.0L - t*(0.125L - t* (251.0L/2880.0L - t*(19.0L/288.0L - t * (19087.0L/362880.0L - t * (751.0L/17280.0L - t * (1070017.0L/29030400.0L - t * (2857.0L/89600.0L)))))))));
        if ((r >= m+s)){ // if (r - m - t >= 0.0l)
            t = exp(r-m-s);
            lead = r - mterm -log_e_x_m1(t);
        }else{
            lead = s;
            s = exp(r-m-s);
            lead += -mterminv - s * (0.5 + s * (1.0/24.0 - s * s *(1.0/2880.0 + s * (1.0/840 - s * (217.0/508032 -s * (635229.0/5080320 - s* (5039857.0/29030400)))))));
        }
    }
    mterm = log_e_x_p1(m);

    fout =0.0;

    if (r >= 18.375f){
        s = exp(-r);
        t = s + 1.0;
        tmptmp = lead;
        for(uint32_t i=1;i<x;i++, t+= s) {
            fout += exp(tmptmp -lngamma(1.0+i));
            tmptmp += r + log(t) - mterm;
        }
        fout += 0.5 * exp(tmptmp -lngamma(1.0+x));
    }else{
        s = exp(r);
        t = 1.0 + s;// F(r) = sum  exp(r) *k / (k + exp(r))*(k + exp(r))
        tmptmp = lead;
        for(uint32_t i=1;i<x;i++, t+= 1.0) {
            fout += exp(tmptmp -lngamma(1.0+i));
            tmptmp += log(t) - mterm;
        }
        fout += 0.5 * exp(tmptmp -lngamma(1.0+x) );
    }

return fout;
}

void deflatedNB_Zscore_d2(double* fout6, const uint32_t x, const double r, const double m){
    double s,t;
    double mterm = log_e_x_p1(m);
    double mterminv = log_e_x_p1(-m);
    double dbuf6[6];

    if (m <= 6.0){
        if (r >= -12.0f){
            dbuf6[0] = r - mterm -log_e_x_m1(mterminv *exp(r));
            s = exp(r) /(1.0 - exp(-mterminv *exp(r)));
            dbuf6[1] = 1.0 - mterminv * s;
            dbuf6[2] = (s-exp(m)) / (exp(m)+1);
            dbuf6[3] = s * mterminv *(s * exp(-mterminv *exp(r)) * mterminv - 1.0);
            t = exp(-mterminv *exp(r));
            dbuf6[4] = exp(m) * (-exp(r)-1.0 - t * (t - exp(2.0*r-m) - 2.0 - exp(r) )) / ((t-1)*(t-1)*(exp(m)+1)*(exp(m)+1));
            dbuf6[5] = (s / (exp(m)+1.0)) * (1.0 - mterminv * exp(r) / ( (-1.0 + exp(mterminv *exp(r)))));
            //dbuf6[5] = (s / (exp(m)+1)) * s * (1.0 + mterminv / (-1.0 + exp((log(exp(m)+1.0)-m) *exp(r))));
        }else{
            t = mterminv;
            s = t *exp(r);
            dbuf6[0] = -log(mterminv) - mterm - 0.5 * s;
            dbuf6[1] = -s * (0.5 + s * (2.0/24.0 - s * s *(4.0/2880.0 + s * (5.0/840 - s * (1302.0/508032 -s * (635229.0*7.0/5080320 - s* (8.0*5039857.0/29030400)) )))));
            dbuf6[2] = (1.0 - exp(m) * t) / ((exp(m)+1)*t);
            dbuf6[3] = -s * (0.5 - s * (4.0/24.0 - s * s *(16.0/2880.0 - s * (25.0/840 + s * (6.0*1302.0/508032 +s * (635229.0*49.0/5080320 + s* (5039857.0*64.0/29030400)) )))));
            dbuf6[4] = (1.0 - exp(m) * t *(t + 1.0)) / (t * t * (exp(m) +1) * (exp(m) +1));
            dbuf6[5] = -dbuf6[3] / (t * (exp(m)+1));
        }
    }else{
        t = exp(-m);
        s = t * (0.5L - t * (5.0L/24.0L - t*(0.125L - t* (251.0L/2880.0L - t*(19.0L/288.0L - t * (19087.0L/362880.0L - t * (751.0L/17280.0L - t * (1070017.0L/29030400.0L - t * (2857.0L/89600.0L)))))))));
        dbuf6[2] = t * (0.5 - t * (5.0/12.0 - t*(0.375 - t* (251.0/720.0 - t*(95.0/288.0 - t * (19087.0/60480.0 - t * (7.0*751.0/17280.0 - t * (8.0*1070017.0/29030400.0 - t * (9.0*2857.0/89600.0)))))))));
        if ((r >= m+s)){ // if (r - m - t >= 0.0l)
            dbuf6[4] = t * (0.5 - t * (10.0/12.0 - t*(0.375*3.0 - t* (251.0*4.0/720.0 - t*(95.0*5.0/288.0 - t * (36.0*19087.0/362880.0 - t * (49.0*751.0/17280.0 - t * (64.0*1070017.0/29030400.0 - t * (81.0*2857.0/89600.0)))))))));
            t = exp(r-m-s);
            dbuf6[0] = r - mterm -log_e_x_m1(t);
            dbuf6[1] = -t / (1.0 - exp(-t));
            //dbuf6[4] = 1.0 + 1.0/(exp(-m)+1) + t * (dbuf6[4] + (dbuf6[2]-1))/ (1.0 - exp(-t)) - exp(-t) * t * t * (dbuf6[2]-1)/ ((1.0 - exp(-t))*(1.0 - exp(-t)));
            //-1.0 / (1.0 + exp(m)) + dbuf6[2] - t * (0.5 + t * (4.0/24.0 - t * t *(16.0/2880.0 + t * (25.0/840 - t * (6.0*1302.0/508032 -t * (635229.0*49.0/5080320 - t* (64.0*5039857.0/29030400)) )))));
             dbuf6[3] = dbuf6[1] + exp(-t) * t * t / ((1.0 - exp(-t))*(1.0 - exp(-t)));
            dbuf6[5] = (dbuf6[1] + exp(-t) * t * t / ((1.0 - exp(-t))*(1.0 - exp(-t)))) * (-1.0 - dbuf6[2]);
            dbuf6[4] = -exp(-m) / ((exp(-m)+1)*(exp(-m)+1)) - dbuf6[4] * dbuf6[1] + (fout6[2]-1.0) * fout6[5];
            dbuf6[2] = -1.0 / (exp(-m)+1) + (dbuf6[2]-1.0) * dbuf6[1];
            dbuf6[1] += 1.0;
               //        -exp(r-m-t) / (1.0 - exp(-exp(r-m-t))) + exp(-exp(r-m-t)) * exp(2.0 * (r-m-t)) / ((1.0 - exp(-exp(r-m-t)))*(1.0 - exp(-exp(r-m-t)))); // TOCHECK
        }else{
            dbuf6[0] = s;
            s = exp(r-m-s);
            dbuf6[0] += -mterminv - s * (0.5 + s * (1.0/24.0 - s * s *(1.0/2880.0 + s * (1.0/840 - s * (217.0/508032 -s * (635229.0/5080320 - s* (5039857.0/29030400)))))));
            dbuf6[1] = -s * (0.5 + s * (2.0/24.0 - s * s *(4.0/2880.0 + s * (5.0/840 - s * (1302.0/508032 -s * (635229.0*7/5080320 - s* (8*5039857.0/29030400)) )))));
            dbuf6[4] = -1.0 / (1.0 + exp(m)) + dbuf6[2] - s * (0.5 + s * (4.0/24.0 - s * s *(16.0/2880.0 + s * (25.0/840 - s * (6.0*1302.0/508032 -s * (635229.0*49.0/5080320 - s* (64.0*5039857.0/29030400)) )))));
            //1.0 + 1.0/(exp(-m)+1) + s * (dbuf6[4] + (dbuf6[2]-1))/ (1.0 - exp(-s)) - exp(-s) * s * s * (dbuf6[2]-1)/ ((1.0 - exp(-s))*(1.0 - exp(-s)));
            dbuf6[3] = -s * (0.5 + s * (4.0/24.0 - s * s *(16.0/2880.0 + s * (25.0/840 - s * (6.0*1302.0/508032 -s * (635229.0*49.0/5080320 - s* (5039857.0*64.0/29030400)) )))));
            dbuf6[5] = -dbuf6[3] * (1.0 + dbuf6[2]);
            dbuf6[2] = 1.0 / (1.0 + exp(m)) - dbuf6[2] - dbuf6[1];
        }
    }


    fout6[0] =0.0;
    fout6[1] =0.0;
    fout6[2] =0.0;
    fout6[3] =0.0;
    fout6[4] =0.0;
    fout6[5] =0.0;

    double tmptmp,lead;
    if (r >= 18.375f){
        s = exp(-r);
        t = s + 1.0;
        tmptmp = lead;
        for(uint32_t i=1;i<x;i++, t+= s) {
            fout6[0] += exp(tmptmp -lngamma(1.0+i));
            tmptmp += r + log(t) - mterm;
        }
        fout6[0] += 0.5 * exp(tmptmp -lngamma(1.0+x));
    }else{
        s = exp(r);
        t = 1.0 + s;// F(r) = sum  exp(r) *k / (k + exp(r))*(k + exp(r))
        tmptmp = lead;
        for(uint32_t i=1;i<x;i++, t+= 1.0) {
            fout6[0] += exp(tmptmp -lngamma(1.0+i));
            tmptmp += log(t) - mterm;
        }
        fout6[0] += 0.5 * exp(tmptmp -lngamma(1.0+x) );
    }

    fout6[0] = LogPvalue_to_stdnorm(fout6[0]);
}

uint32_t deflatedNB_sample(const double r, const double m, uint32_t maxval){
    double s,t;
    double lead, fout, tmptmp;
    double mterm = log_e_x_p1(m);
    double mterminv = log_e_x_p1(-m);
    if (m <= 3.0){
        if (r >= -12.0f){
            lead = r - mterm -log_e_x_m1(mterminv *exp(r));
        }else{
            t = mterminv;
            s = t *exp(r);
            lead = -log(mterminv) -log(exp(m)+1.0) - 0.5 * s;
        }
    }else{
        t = exp(-m);
        s = t * (0.5L - t * (5.0L/24.0L - t*(0.125L - t* (251.0L/2880.0L - t*(19.0L/288.0L - t * (19087.0L/362880.0L - t * (751.0L/17280.0L - t * (1070017.0L/29030400.0L - t * (2857.0L/89600.0L)))))))));
        if ((r >= m+s)){ // if (r - m - t >= 0.0l)
            t = exp(r-m-s);
            lead = r - mterm -log_e_x_m1(t);
        }else{
            lead = s;
            s = exp(r-m-s);
            lead += -mterminv - s * (0.5 + s * (1.0/24.0 - s * s *(1.0/2880.0 + s * (1.0/840 - s * (217.0/508032 -s * (635229.0/5080320 - s* (5039857.0/29030400)))))));
        }
    }
    mterm = log_e_x_p1(m);

    fout =0.0;
    uint32_t rn; ExOp::toRand(rn);
    double pval = pow(0.5,32.0) * rn;
    uint32_t x;
    if (r >= 18.375f){
        s = exp(-r);
        t = s + 1.0;
        tmptmp = lead;
        for(x=1;((pval>0.0)&&(x<maxval));x++, t+= s) {
            pval -= exp(tmptmp -lngamma(1.0+x));
            tmptmp += r + log(t) - mterm;
        }
    }else{
        s = exp(r);
        t = 1.0 + s;// F(r) = sum  exp(r) *k / (k + exp(r))*(k + exp(r))
        tmptmp = lead;
        for(x=1;((pval>0.0)&&(x<maxval));x++, t+= 1.0) {
            pval -= exp(tmptmp -lngamma(1.0+x));
            tmptmp += log(t) - mterm;
        }
    }
return x;}
void deflatedNB_getMeanAndVar(double* fout2, const double r, const double m){
    uint32_t x;
    double buffer[6];
   /* WeightElem<double,2u> hehe;
    hehe.toZero();
    for(x=1 ; x< 25;x++){
        deflatedNBterm_wr_f_d2(buffer,x,r,m);
        hehe += WeightElem<double,2u>((double)x, exp(buffer[0] - lngamma(1.0 + x)) );
    }


    if (hehe.w[0] > 0.9f) {
        fout2[0] = hehe.getMean();
        fout2[1] = hehe.getVar_biaised();
    }else{*/
        double p = 1.0 / (1.0 + exp(m));
        double rr = exp(r);
        double denum;
        if (m <= 3.0){
            if (r < -12.0f) rr = exp(-12.0);
            denum = -rr * (log(1.0 + exp(m)) - m);
            denum = (1.0 - exp(denum  )) / rr;
            fout2[0] = 1.0 / (denum * exp(m));
            fout2[1] = 1.0 * ((1.0 + rr) + exp(m))  / (exp(2*m) * denum) - fout2[0]* fout2[0];
        }else{
            /*if (r >= 12.0f){
                denum = -rr * exp(-m) * log_e_x_1plus_overe_x(m);
                denum = log_e_x_1plus_overe_x(-m) * e_x_1minus_overx(denum);
                fout2[0] = 1.0 / (denum);
                fout2[1] = (exp(-m) * (1.0 + rr) + 1.0) / (denum) - fout2[0]* fout2[0];
            }else{*/
                denum = -rr * log(1.0 + exp(-m));
                denum = log_e_x_1plus_overe_x(-m) * e_x_1minus_overx(denum);
                fout2[0] = 1.0 / (denum);
                fout2[1] = (exp(-m) * (1.0 + rr) + 1.0) / (denum) - fout2[0]* fout2[0];
            //}
        }
    //}
}
void deflatedNB_getMeanAndVarSilly(double* fout2, const double r, const double m){
    uint32_t x;
    double buffer[6];
    WeightElem<double,2u> hehe;
    hehe.toZero();
    for(x=1 ; x< 25;x++){
        deflatedNBterm_wr_f_d2(buffer,x,r,m);
        hehe += WeightElem<double,2u>((double)x, exp(buffer[0] - lngamma(1.0 + x)) );
    }
    deflatedNB_getMeanAndVar(fout2, r,m);
    fout2[0] -= hehe.getMean();
    fout2[1] -= hehe.getVar_biaised();
}


double deflatedNB_pvalueSilly(const uint32_t v, const double r, const double m){
    uint32_t x;
    double buffer[6];
    double hehe = 0.0;
    for(x=1 ; x< v;x++){
        deflatedNBterm_wr_f_d2(buffer,x,r,m);
        hehe += exp(buffer[0] - lngamma(1.0 + x));
    }
    deflatedNBterm_wr_f_d2(buffer,v,r,m);
    hehe += 0.5 * exp(buffer[0] - lngamma(1.0 + v));
return hehe;}

LFH_GOLD double LogProb_NBdistrib_exppara(uint32_t x, double log_r, double logit_p){double r =exp(log_r); return lngammachoose(x, r - 1.0 +  x) -log(1.0 + exp(logit_p)) * r - log(1.0 + exp(-logit_p)) * x;}

#ifdef GNU_SCIENTIFIC_LIBRARY
#include "GSLfunc.hpp"

static double log_beta_logit_cont_frac(const double a, const double b, const double x){
	const double cutoff = 2.0 * ExCo<double>::delta();  /* control the zero cutoff */
	unsigned int iter_count = 0;
	double cf;

	if (x < -18.0f){
        return (a+b) * exp(x) / (a+1.0);
	}

    double lx = 1.0 / (1.0 + exp(-x));


	/* standard initialization for continued fraction */
	double num_term = 1.0;
	double den_term = 1.0 - (a+b)*lx/(a+1.0);
	if (fabs(den_term) < cutoff) den_term = cutoff;
	den_term = 1.0/den_term;
	cf = log(den_term);

	while(iter_count < 256u) {
		const int k  = iter_count + 1;
		double coeff = k*(b-k)*lx/(((a-1.0)+2*k)*(a+2*k));
		double delta_frac;

		/* first step */
		den_term = 1.0 + coeff*den_term;
		num_term = 1.0 + coeff/num_term;
		if(fabs(den_term) < cutoff) den_term = cutoff;
		if(fabs(num_term) < cutoff) num_term = cutoff;
		den_term  = 1.0/den_term;

		delta_frac = den_term * num_term;
		cf += log(delta_frac);

		coeff = -(a+k)*(a+b+k)*lx/((a+2*k)*(a+2*k+1.0));

		/* second step */
		den_term = 1.0 + coeff*den_term;
		num_term = 1.0 + coeff/num_term;
		if(fabs(den_term) < cutoff) den_term = cutoff;
		if(fabs(num_term) < cutoff) num_term = cutoff;
		den_term = 1.0/den_term;

		delta_frac = den_term*num_term;
		cf += log(delta_frac);

		if(fabs(delta_frac-1.0) < 2.0*GSL_DBL_EPSILON) break;

		++iter_count;
	}
return cf;}

double lfh_sf_log_beta_logit_inc_e(const double a, const double b, const double x){
	if (a <= 0.0 || b <= 0.0) return(NAN);
	else if (x == -INFINITY) return(-INFINITY);
	else if (x == INFINITY) return 0;
    const double ln_pre_val = -lnbeta(a,b) - a * log_e_x_p1(-x) - b * log_e_x_p1(x);
    if (exp(x)*(b + 1.0) < (a+1.0)) return ln_pre_val + log_beta_logit_cont_frac(a, b, x) - log(a);
    /* Apply continued fraction after hypergeometric transformation. */
    const double term = ln_pre_val + log_beta_logit_cont_frac(b, a, -x) - log(b);
    return log_1_me_x(term);
}

static double log_beta_exp_logit_cont_frac(const double a, const double b, const double x){ //
	const double cutoff = 2.0 * ExCo<double>::delta();  /* control the zero cutoff */
	unsigned int iter_count = 0;
	double cf;

    double ea = x +log_e_x_p1(b-a) - log_e_x_p1(-a);

	if (ea < -18.0f) return exp(ea) - a;
	ea = exp(a);
	double eb = exp(b);

    double lx = 1.0 / (1.0 + exp(-x));

	/* standard initialization for continued fraction */
	double num_term = 1.0;
	double den_term = 1.0 - (ea+eb)*lx/(ea+1.0);
	if (fabs(den_term) < cutoff) den_term = cutoff;
	den_term = 1.0/den_term;
	cf = log(den_term);

	while(iter_count < 256u) {
		const int k  = iter_count + 1;
		double coeff = k*(eb-k)*lx/(((ea-1.0)+2*k)*(ea+2*k));
		double delta_frac;

		/* first step */
		den_term = 1.0 + coeff*den_term;
		num_term = 1.0 + coeff/num_term;
		if(fabs(den_term) < cutoff) den_term = cutoff;
		if(fabs(num_term) < cutoff) num_term = cutoff;
		den_term  = 1.0/den_term;

		delta_frac = den_term * num_term;
		cf += log(delta_frac);

		coeff = -(ea+k)*(ea+eb+k)*lx/((ea+2*k)*(ea+2*k+1.0));

		/* second step */
		den_term = 1.0 + coeff*den_term;
		num_term = 1.0 + coeff/num_term;
		if(fabs(den_term) < cutoff) den_term = cutoff;
		if(fabs(num_term) < cutoff) num_term = cutoff;
		den_term = 1.0/den_term;

		delta_frac = den_term*num_term;
		cf += log(delta_frac);

		if(fabs(delta_frac-1.0) < 2.0*GSL_DBL_EPSILON) break;

		++iter_count;
	}
return cf - a;}

double lfh_sf_log_beta_exp_logit_inc_e(const double a, const double b, const double x){
    if (x == -INFINITY) return(-INFINITY);
	else if (x == INFINITY) return 0;
	double ea = exp(a);
	double eb = exp(b);
    const double ln_pre_val = -lnbetaexp(a,b) - ea * log_e_x_p1(-x) - eb * log_e_x_p1(x);
    if (exp(x)*(eb+1.0) < (ea+1.0)) return ln_pre_val + log_beta_exp_logit_cont_frac(a, b, x);
    else return log_1_me_x(ln_pre_val + log_beta_exp_logit_cont_frac(b, a, -x));
}


// logit_transformed input "x", needs to use logistic function to recover




double deflatedNB_pvalue(const uint32_t x, const double r, const double m){
    double s,t;
    double lead, fout, tmptmp;
    //double mterm = log_e_x_p1(m);
    double mterminv = log_e_x_p1(-m);
    /*if (m <= 3.0){
        if (r >= -12.0f){
            lead = r - mterm -log_e_x_m1(mterminv *exp(r));
        }else{
            t = mterminv;
            s = t *exp(r);
            lead = -log(mterminv) -log(exp(m)+1.0) - 0.5 * s;
        }
    }else{
        t = exp(-m);
        s = t * (0.5L - t * (5.0L/24.0L - t*(0.125L - t* (251.0L/2880.0L - t*(19.0L/288.0L - t * (19087.0L/362880.0L - t * (751.0L/17280.0L - t * (1070017.0L/29030400.0L - t * (2857.0L/89600.0L)))))))));
        if ((r >= m+s)){ // if (r - m - t >= 0.0l)
            t = exp(r-m-s);
            lead = r - mterm -log_e_x_m1(t);
        }else{
            lead = s;
            s = exp(r-m-s);
            lead += -mterminv - s * (0.5 + s * (1.0/24.0 - s * s *(1.0/2880.0 + s * (1.0/840 - s * (217.0/508032 -s * (635229.0/5080320 - s* (5039857.0/29030400)))))));
        }
    }*/

    double zerop = log(0.5) - log_1_me_x(-mterminv * exp(r));
    double tmp = lfh_sf_log_beta_exp_logit_inc_e(log((double)(x)),r,-m);
    return 1.0 - exp(tmp+ zerop + log_e_x_p1(lfh_sf_log_beta_exp_logit_inc_e(log((double)(x+1)),r,-m)-tmp));
}

double deflatedNB_logpvalue(const uint32_t x, const double r, const double m){
    double s,t;
    double lead, fout, tmptmp;
    //double mterm = log_e_x_p1(m);
    double mterminv = log_e_x_p1(-m);
    /*if (m <= 3.0){
        if (r >= -12.0f){
            lead = r - mterm -log_e_x_m1(mterminv *exp(r));
        }else{
            t = mterminv;
            s = t *exp(r);
            lead = -log(mterminv) -log(exp(m)+1.0) - 0.5 * s;
        }
    }else{
        t = exp(-m);
        s = t * (0.5L - t * (5.0L/24.0L - t*(0.125L - t* (251.0L/2880.0L - t*(19.0L/288.0L - t * (19087.0L/362880.0L - t * (751.0L/17280.0L - t * (1070017.0L/29030400.0L - t * (2857.0L/89600.0L)))))))));
        if ((r >= m+s)){ // if (r - m - t >= 0.0l)
            t = exp(r-m-s);
            lead = r - mterm -log_e_x_m1(t);
        }else{
            lead = s;
            s = exp(r-m-s);
            lead += -mterminv - s * (0.5 + s * (1.0/24.0 - s * s *(1.0/2880.0 + s * (1.0/840 - s * (217.0/508032 -s * (635229.0/5080320 - s* (5039857.0/29030400)))))));
        }
    }*/

    double zerop = log(0.5) - log_1_me_x(-mterminv * exp(r));
    double tmp = lfh_sf_log_beta_exp_logit_inc_e(log((double)(x)),r,-m);
    return log_1_me_x(tmp+ zerop + log_e_x_p1(lfh_sf_log_beta_exp_logit_inc_e(log((double)(x+1)),r,-m)-tmp));
}

double NB_pvalue(const uint32_t x, const double r, const double m){
    double tmp = lfh_sf_log_beta_exp_logit_inc_e(log((double)(x)),r,-m);
    return 1.0 - exp(tmp+ log(0.5) + log_e_x_p1(lfh_sf_log_beta_exp_logit_inc_e(log((double)(x+1)),r,-m)-tmp));
}

double NB_logpvalue(const uint32_t x, const double r, const double m){
    double tmp = lfh_sf_log_beta_exp_logit_inc_e(log((double)(x)),r,-m);
    return log_1_me_x(tmp+ log(0.5) + log_e_x_p1(lfh_sf_log_beta_exp_logit_inc_e(log((double)(x+1)),r,-m)-tmp));
}


double NB_logitpvalue(const uint32_t x, const double r, const double m, const double dropout_term){
    double pval = deflatedNB_pvalue(x, r, m);
    double logrp0 = exp(r) * log_e_x_p1(m); // probability of zero;
    double posterior_p0 = exp(logrp0) * log_e_x_p1(dropout_term);
    double trval = pval * (1.0 - posterior_p0) + posterior_p0;
    return log(trval) - log(1.0 - trval);
}


LFH_GOLD double d_lngamma_dx(double x){gsl_sf_result_struct res;gsl_sf_psi_e(x,&res); return res.val;}
LFH_GOLD double d2_lngamma_dx2(double x){gsl_sf_result_struct res;gsl_sf_psi_1_e(x,&res); return res.val;}
/*LFH_GOLD double d_lngamma_dx_intdiff(double x, int k){
    if (abs(k) <= 1){
        switch(k){
        case -1: return 1.0 / (1.0 - x);
        case 0: return 0.0;
        case 1: return 1.0 / x;
        }
    }else return d_lngamma_dx(x + k) - d_lngamma_dx(x);
}*/
LFH_GOLD double Pvalue_to_stdnorm(double x){return sliwa_gsl_cdf_ugaussian_Pinv(x);}
LFH_GOLD double LogPvalue_to_stdnorm(double x){return handfield_gsl_cdf_ugaussian_logPinv(x);}

LFH_GOLD double Pvalue_chisquarre_Ptail(double x, double free){return incgamma_frac(free* 0.5, x * 0.5);}
LFH_GOLD double Pvalue_chisquarre_Ntail(double x, double free){return 1.0 - incgamma_frac(free* 0.5, x* 0.5);}
LFH_GOLD double Pvalue_Gamma_Ptail(double x, double k, double theta){return incgamma_frac(k, x / theta);}
LFH_GOLD double Pvalue_Gamma_Ntail(double x, double k, double theta){return 1.0 - incgamma_frac(k, x / theta);}
LFH_GOLD double Pvalue_GammaRate_Ptail(double x, double alpha, double beta){return incgamma_frac(alpha, x* beta);}
LFH_GOLD double Pvalue_GammaRate_Ntail(double x, double alpha, double beta){return 1.0 - incgamma_frac(alpha, x* beta);}
LFH_GOLD double Pvalue_SumLogPvalues_Ptail(double sum, int nbpvals){ return incgamma_frac((double) nbpvals , -sum); }

LFH_GOLD double Pvalue_Beta_Ntail(double x, double a, double b){gsl_sf_result_struct res;gsl_sf_beta_inc_e(b,a,1.0-x,&res);	return res.val;}
LFH_GOLD double Pvalue_Beta_Ptail(double x, double a, double b){gsl_sf_result_struct res;gsl_sf_beta_inc_e(a,b,x,&res);	return res.val;}
LFH_GOLD double Pvalue_Fdistrib_Ntail(double x, double a, double b){gsl_sf_result_struct res;gsl_sf_beta_inc_e(0.5*a,0.5*b,a*x/(a*x+b),&res); return res.val;}
LFH_GOLD double Pvalue_Fdistrib_Ptail(double x, double a, double b){gsl_sf_result_struct res;gsl_sf_beta_inc_e(0.5*b,0.5*a,b/(a*x+b),&res);	return res.val;}

LFH_GOLD double Pvalue_Tdistrib_Ptail(double x, double v){return gsl_cdf_tdist_Q(x,v);}
LFH_GOLD double Pvalue_Tdistrib_Ntail(double x, double v){return Pvalue_Tdistrib_Ptail(-x, v);}
LFH_GOLD double Pvalue_to_Tdistrib(double x, double v){return gsl_cdf_tdist_Pinv(x,v);}
LFH_GOLD double logPvalue_to_Tdistrib(double x, double v){return handfield_gsl_cdf_tdist_logPinv(x,v);}



// log are all natural logarithms

LFH_GOLD double LogPvalue_chisquarre_Ptail(double x, double free){return incgamma_log_frac(free* 0.5, x * 0.5);}
LFH_GOLD double LogPvalue_chisquarre_Ntail(double x, double free){return incgamma_log_1minus_frac(free* 0.5, x* 0.5);}
LFH_GOLD double LogPvalue_Gamma_Ptail(double x, double k, double theta){return incgamma_log_frac(k, x / theta);}
LFH_GOLD double LogPvalue_Gamma_Ntail(double x, double k, double theta){return incgamma_log_1minus_frac(k, x / theta);}
LFH_GOLD double LogPvalue_GammaRate_Ptail(double x, double alpha, double beta){return incgamma_log_frac(alpha, x* beta);}
LFH_GOLD double LogPvalue_GammaRate_Ntail(double x, double alpha, double beta){return incgamma_log_1minus_frac(alpha, x* beta);}
LFH_GOLD double LogPvalue_SumLogPvalues_Ptail(double sum, int nbpvals){ return incgamma_log_frac((double) nbpvals ,-sum); }

LFH_GOLD double LogPvalue_Beta_Ntail(double x, double a, double b){gsl_sf_result_struct res;gsl_sf_log_beta_inc_e(b,a,1.0-x,&res);	return res.val;}
LFH_GOLD double LogPvalue_Beta_Ptail(double x, double a, double b){gsl_sf_result_struct res;gsl_sf_log_beta_inc_e(a,b,x,&res);	return res.val;}
LFH_GOLD double LogPvalue_Fdistrib_Ntail(double x, double a, double b){gsl_sf_result_struct res;gsl_sf_log_beta_inc_e(0.5*a,0.5*b,a*x/(a*x+b),&res); return res.val;}
LFH_GOLD double LogPvalue_Fdistrib_Ptail(double x, double a, double b){gsl_sf_result_struct res;gsl_sf_log_beta_inc_e(0.5*b,0.5*a,b/(a*x+b),&res); return res.val;}
	/*
double Gcummul_HotelingTTest(double stat, unsigned int nb_dims, unsigned int nbsample_A, unsigned int nbsample_b){
	// stat = [(NaNb) * (Na+Nb-D-1) / (D*(Na+Nb))] *  ( ((sum_allA X)/Na) - ((sum_allB X)/Nb))^T ( sum_allA&B XX^T )^-1 (((sum_allA X)/Na) - ((sum_allB X)/Nb))
}*/

LFH_GOLD	double Pvalue_NBdistrib_GT(uint32_t x, double r, double p){return exp(LogPvalue_NBdistrib_GT(x,r,p));}
LFH_GOLD	double Pvalue_NBdistrib_LE(uint32_t x, double r, double p){return exp(LogPvalue_NBdistrib_LE(x,r,p));}
LFH_GOLD	double Pvalue_NBdistrib_GE(uint32_t x, double r, double p){return (x == 0) ? 1.0 : exp(LogPvalue_NBdistrib_GE(x,r,p));}
LFH_GOLD	double Pvalue_NBdistrib_LT(uint32_t x, double r, double p){return (x == 0) ? 0.0 : exp(LogPvalue_NBdistrib_LT(x,r,p));}
LFH_GOLD	double Pvalue_NBdistrib_GH(uint32_t x, double r, double p){return exp(LogPvalue_NBdistrib_GH(x,r,p));}
LFH_GOLD	double Pvalue_NBdistrib_LH(uint32_t x, double r, double p){return exp(LogPvalue_NBdistrib_LH(x,r,p));}
LFH_GOLD	double LogPvalue_NBdistrib_GT(uint32_t x, double r, double p){gsl_sf_result_struct res;gsl_sf_log_beta_inc_e(1.0 + x,r, p ,&res); return res.val;}
LFH_GOLD	double LogPvalue_NBdistrib_LE(uint32_t x, double r, double p){gsl_sf_result_struct res;gsl_sf_log_beta_inc_e(r,1.0 + x, 1.0 - p ,&res);return res.val;}
LFH_GOLD	double LogPvalue_NBdistrib_GE(uint32_t x, double r, double p){gsl_sf_result_struct res;gsl_sf_log_beta_inc_e((double)x,r, p ,&res); return res.val;}
LFH_GOLD	double LogPvalue_NBdistrib_LT(uint32_t x, double r, double p){gsl_sf_result_struct res;gsl_sf_log_beta_inc_e(r,(double)x, 1.0 - p ,&res);return res.val;}
LFH_GOLD	double LogPvalue_NBdistrib_GH(uint32_t x, double r, double p){gsl_sf_result_struct res;gsl_sf_log_beta_inc_e(0.5 + x,r, p ,&res); return res.val;}
LFH_GOLD	double LogPvalue_NBdistrib_LH(uint32_t x, double r, double p){gsl_sf_result_struct res;gsl_sf_log_beta_inc_e(r,0.5 + x, 1.0 - p ,&res);return res.val;}

//LFH_GOLD	double LogProb_NBdistrib_exppara(uint32_t x, double log_r, double logit_p){double r = exp(log_r); return return log(1.0 + exp(logit_p)) * r + log(1.0 + exp(-logit_p)) * x + lngammachoose(r - 1.0 +  x , x);}
LFH_GOLD	double LogPvalue_NBdistrib_exppara_GT(uint32_t x, double log_r, double logit_p){return lfh_gsl_sf_log_beta_logit_inc_e(1.0 + x,exp(log_r), logit_p);}
LFH_GOLD	double LogPvalue_NBdistrib_exppara_LE(uint32_t x, double log_r, double logit_p){return lfh_gsl_sf_log_beta_logit_inc_e(exp(log_r),1.0 + x,-logit_p);}
LFH_GOLD	double LogPvalue_NBdistrib_exppara_GE(uint32_t x, double log_r, double logit_p){return lfh_gsl_sf_log_beta_logit_inc_e(x,exp(log_r), logit_p);}
LFH_GOLD	double LogPvalue_NBdistrib_exppara_LT(uint32_t x, double log_r, double logit_p){return lfh_gsl_sf_log_beta_logit_inc_e(exp(log_r),x,-logit_p);}

LFH_GOLD	double LogPvalue_NBdistrib_exppara_GH(uint32_t x, double log_r, double logit_p){return lfh_gsl_sf_log_beta_logit_inc_e(0.5 + x,exp(log_r), logit_p);}
LFH_GOLD	double LogPvalue_NBdistrib_exppara_LH(uint32_t x, double log_r, double logit_p){return lfh_gsl_sf_log_beta_logit_inc_e(exp(log_r),0.5 + x,-logit_p);}
LFH_GOLD	double LogPvalue_NBdistrib_exppara_GH_exact(uint32_t x, double log_r, double logit_p){double res[2]; res[0] = lfh_gsl_sf_log_beta_logit_inc_e(x,exp(log_r), logit_p); res[1] = lfh_gsl_sf_log_beta_logit_inc_e(1.0 + x,exp(log_r), logit_p); return res[0] + (log_e_x_p1(res[1] - res[0])- M_LN2);}
LFH_GOLD	double LogPvalue_NBdistrib_exppara_LH_exact(uint32_t x, double log_r, double logit_p){double res[2]; res[0] = lfh_gsl_sf_log_beta_logit_inc_e(exp(log_r),1.0 + x,-logit_p); res[1] = lfh_gsl_sf_log_beta_logit_inc_e(exp(log_r),x,-logit_p); return res[0] + (log_e_x_p1(res[1] - res[0])- M_LN2);}






//double mean_of_clamped_gamma2(double k, double theta, double max){ // mean of gamma distributed random variable, conditionned that it cannot be higher than *max*
//	return (max <=0) ? max : theta * k * incgamma_1minus_frac(k+1,max/theta) / incgamma_1minus_frac(k,max/theta);
//	}

double mean_of_clamped_gamma(double k, double theta, double max){ // mean of gamma distributed random variable, conditionned that it cannot be higher than *max*
//	return (max <=0) ? max : theta * gamma(k+1)*incgamma_1minus_frac(k+1,max/theta) / gamma(k)*incgamma_1minus_frac(k,max/theta);
//	return (max <=0) ? max : theta * (k * gamma(k)*incgamma_1minus_frac(k,max/theta) - pow(max/theta,k)*exp(-max/theta)) / gamma(k)*incgamma_1minus_frac(k,max/theta);
	if (max <=0.0f) return max;
	gsl_sf_result gres;
	double fout =-log(k) -(max/theta) -lngamma(k)+ k * log(max/theta)-incgamma_log_1minus_frac(k,max/theta);
	if (fout >= 0.0f) return 0.0f;
	if (!ExOp::isValid(fout)) return max;
	gsl_sf_expm1_e( fout ,&gres);
	gres.val = -gres.val * theta * k;
	return (ExOp::isValid(gres.val) && (gres.val<= max)) ?( (gres.val < 0) ? 0.0f : gres.val) : max;
    }
double mean_of_clamped_gamma_rate(double k, double beta, double max){ // mean of gamma distributed random variable, conditionned that it cannot be higher than *max*
	//	return (max <=0) ? max : theta * gamma(k+1)*incgamma_1minus_frac(k+1,max/theta) / gamma(k)*incgamma_1minus_frac(k,max/theta);
	//	return (max <=0) ? max : theta * (k * gamma(k)*incgamma_1minus_frac(k,max/theta) - pow(max/theta,k)*exp(-max/theta)) / gamma(k)*incgamma_1minus_frac(k,max/theta);
	if (max <=0.0f) return max;
	gsl_sf_result gres;
	double fout =-log(k) -(max*beta) -lngamma(k)+ k * log(max*beta)-incgamma_log_1minus_frac(k,max*beta);
	if (fout >= 0.0f) return 0.0f;
	if (!ExOp::isValid(fout)) return max;
	gsl_sf_expm1_e( fout ,&gres);
	gsl_sf_expm1_e( fout ,&gres);
	gres.val = -gres.val *  k / beta;
	return (ExOp::isValid(gres.val) && (gres.val<= max)) ?( (gres.val < 0) ? 0.0f : gres.val) : max;
   }




#endif

double log_beta_cont_frac(const double a, const double b, const double x){
	const unsigned int max_iter = 512;        // control iterations
	const double cutoff = 2.0 * ExCo<double>::delta();  // control the zero cutoff
	unsigned int iter_count = 0;
	double cf;

	// standard initialization for continued fraction
	double num_term = 1.0;
	double den_term = 1.0 - (a+b)*x/(a+1.0);
	if (fabs(den_term) < cutoff) den_term = cutoff;
	den_term = 1.0/den_term;
	cf = log(den_term);
    double delta;
	while(iter_count < max_iter) {
		const int k  = iter_count + 1;
		double coeff = k*(b-k)*x/(((a-1.0)+2*k)*(a+2*k));
		// first step
		den_term = 1.0 + coeff*den_term;
		num_term = 1.0 + coeff/num_term;
		delta =  (fabs(num_term) < cutoff) ? log(cutoff) : log(num_term);
		delta -= (fabs(den_term) < cutoff) ? log(cutoff) : log(den_term);
		den_term  = 1.0/den_term;

		cf += delta;

		coeff = -(a+k)*(a+b+k)*x/((a+2*k)*(a+2*k+1.0));

		// second step
		den_term = 1.0 + coeff*den_term;
		num_term = 1.0 + coeff/num_term;
		delta = (fabs(num_term) < cutoff) ? log(cutoff) : log(num_term);
		delta -= (fabs(den_term) < cutoff) ? log(cutoff) : log(den_term);
        den_term  = 1.0/den_term;

		cf += delta;

		if(fabs(delta) < 2.0*GSL_DBL_EPSILON) break;

		++iter_count;
	}


	if(iter_count >= max_iter){
        printf("Warning, could not converge in time, %e %e %e as input\n", a, b, x);
	}
	return cf;
}
double lfh_gsl_sf_log_beta_logit_inc_e(const double a, const double b, const double logit_x){
	if ((a <= 0.0) || (b <= 0.0)) return(-INFINITY);

    gsl_sf_result ln_beta;
    double ln_x = log_e_x_p1(-logit_x);
    double ln_1mx = log_e_x_p1(logit_x);
    gsl_sf_lnbeta_e(a, b, &ln_beta);
    const double ln_pre_val = -ln_beta.val - (a * ln_x + b * ln_1mx);
    double x = 1.0 / (1.0 + exp(-logit_x));

    if(x < (a + 1.0)/(a+b+2.0)) return(log_beta_cont_frac(a, b, x) + ln_pre_val - log(a));
    else return log_1_me_x(log_beta_cont_frac(b, a , 1.0 / (1.0 + exp(logit_x))) - log(b) + ln_pre_val) ;
}



void GPParamSearch_old(const Vector<KeyElem<double, double> >& vec, double &mean, double &sigma_n, double &sigma_s, double &scale){
	CurvatureSearchScope cs; cs.init(0.0000001f, 3);
	double param[3];
	double deriv[3]; double tmp, like;
	Trianglix<double> covar; covar.setSize(vec.getSize());
	unsigned int i,j,k;
	if (vec.getSize() < 3){
		switch(vec.getSize()){
			case 0: mean =0; sigma_s =0;sigma_n =0;scale =0; return;
			case 1: mean =vec[0].d; sigma_s =0;sigma_n =0;scale =0; return;
			case 2: mean =0.5f * (vec[0].d + vec[1].d); sigma_s =0;sigma_n =0.5f * (vec[0].d - vec[1].d) *(vec[0].d - vec[1].d) ;scale =0; return;
		}
	}else{


    // after normalization, noise <= 1 noise = 1 / 1+ exp(param)
	// MarjLikeli =
	// DMarjLikeli/Da = tr(( K^-1yyTK^-1 - K^-1)  * DK/Da)
	Tuple<double> yyy; yyy.setSize(vec.getSize());
	tmp =0.0f;like =0.0f;scale =0.0f;sigma_n=0.0f;
	for(i=0;i<vec.getSize();i++) {yyy[i] = vec[i].d; tmp += vec[i].d; like += vec[i].d * vec[i].d; sigma_n +=vec[i].k;  scale += vec[i].k*vec[i].k;}
	tmp /= vec.getSize(); sigma_n /= vec.getSize();
	like /= vec.getSize(); scale /= vec.getSize();
	scale -= sigma_n; like -= tmp*tmp; sigma_n = like;
	like = pow(like,-0.5f);
	yyy -= tmp; yyy *= like; // centering!
    mean = tmp;

	double bparam[3];
	double blike;
	// check for scale == 0;
	deriv[2] = 0.0f;
		yyy.show();
	param[0] = 0.0f; param[1] = 0.0f;
	do{
	for(k=0,i=0;i<vec.getSize();i++){
		for(j=0;j<i;j++) {covar.data[k++] = exp(param[0]);}
		covar.data[k++] = exp(param[0]) + scale;
	}
		covar.show(); printf("ldet = %e\n",covar.log_determinant());
	like = (covar.Xformed_inner_product_of_inverse(yyy) + covar.log_determinant() + log(2.0f* M_PI)) * -1.0f;
	Trianglix<double> dainverse = covar.mkInverse();
	Trianglix<double> factrig = dainverse.Xformed_outer_product(yyy) - dainverse;
	covar.cell(0,0) -= scale;
	deriv[1] = factrig.cell(0,0);
	for(i=1;i<vec.getSize();i++) {
		deriv[1] += factrig.cell(i,i);
		covar.cell(i,i) -= scale;
	} deriv[1] *= -1.0625f * exp(param[1]) / (1.0f + 2.0f * exp(param[1]) +  exp(2.0f* param[1]));
	deriv[0] = factrig.trace_of_product(covar) * exp(param[0]);

    printf("f[%e,%e] = %e   deriv = {%e,%e,%e}\n", param[0], param[1], like,deriv[0],deriv[1],deriv[2]);
	}while(!cs.updateAscent(like,param,deriv));

	blike = like;
	bparam[0] = param[0];
	bparam[1] = param[1];
	bparam[2] = NAN;

	for(unsigned int test=0;test<100;test++){
		param[0] = sampleGaussian() * 2.0f;
		param[1] = sampleGaussian() * 2.0f;
		param[2] = sampleGaussian() * 2.0f- log( 2.0f*scale);

	do{
        scale = 1.0625f / (1.0f + exp(param[1]));
		for(k=0,i=0;i<vec.getSize();i++){
			for(j=0;j<i;j++) {tmp = vec[i].k - vec[j].k; covar.data[k++] = exp(param[0] - tmp*tmp * exp(param[2]));}
			covar.data[k++] = exp(param[0]) + scale;
		}
		like = (covar.Xformed_inner_product_of_inverse(yyy) + covar.log_determinant() + log(2.0f* M_PI)) * -1.0f;
		Trianglix<double> dainverse = covar.mkInverse();
		Trianglix<double> factrig = dainverse.Xformed_outer_product(yyy) - dainverse;
		covar.cell(0,0) -= scale;
		deriv[1] = factrig.cell(0,0);

		for(i=1;i<vec.getSize();i++) {
			deriv[1] += factrig.cell(i,i);
			covar.cell(i,i) -= scale;
		} deriv[1] *= -1.0625f * exp(param[1]) / (1.0f + 2.0f * exp(param[1]) +  exp(2.0f* param[1]));
		deriv[0] = factrig.trace_of_product(covar) * exp(param[0]);
		for(k=0,i=0;i<vec.getSize();i++){
			for(j=0;j<i;j++) {tmp = vec[i].k - vec[j].k; covar.data[k++] *= -tmp*tmp;}
			covar.data[k++] = 0.0f;
		}
		deriv[2] = factrig.trace_of_product(covar) * exp(param[2]);
        printf("f[%e,%e,%e] = %e   deriv = {%e,%e,%e}\n", param[0], param[1], param[2], like,deriv[0],deriv[1],deriv[2]);
	}while(!cs.updateAscent(like,param,deriv));

		if ((ExOp::isValid(like))&&(like > blike)){
			printf("improved!\n");
		memcpy(bparam, param, sizeof(double)*3);
		blike = like;
	}
	}
	sigma_s = exp(bparam[0]) * sigma_n;
	sigma_n *= 1.0625f / (1.0f + exp(bparam[1]));
	if (bparam[2] == NAN) scale = 0.0f;
	else scale = exp(param[2]);
	}
}

	void GPParamSearch(const Vector<KeyElem<double, double> >& vec, double &mean, double &sigma_n, double &sigma_s, double &scale){
		CurvatureSearchScope cs;
		double param[3];
		double deriv[3]; double tmp, like;

		unsigned int i,j,k;
		if (vec.getSize() < 3){
			switch(vec.getSize()){
				case 0: mean =0; sigma_s =0;sigma_n =0;scale =0; return;
				case 1: mean =vec[0].d; sigma_s =0;sigma_n =0;scale =0; return;
				case 2: mean =0.5f * (vec[0].d + vec[1].d); sigma_s =0;sigma_n =0.5f * (vec[0].d - vec[1].d) *(vec[0].d - vec[1].d) ;scale =0; return;
			}
		}else{


			// after normalization, noise <= 1 noise = 1 / 1+ exp(param)
			// MarjLikeli =
			// DMarjLikeli/Da = tr(( K^-1yyTK^-1 - K^-1)  * DK/Da)
			Tuple<double> yyy; yyy.setSize(vec.getSize());
			Trianglix<double> covar; covar.setSize(yyy.getSize());
			tmp =0.0f;like =0.0f;scale =0.0f;sigma_n=0.0f;
			for(i=0;i<yyy.getSize();i++) {yyy[i] = vec[i].d; tmp += vec[i].d; like += vec[i].d * vec[i].d; sigma_n +=vec[i].k;  scale += vec[i].k*vec[i].k;}
			tmp /= yyy.getSize(); sigma_n /= yyy.getSize();
			like /= yyy.getSize(); scale /= yyy.getSize();
			scale -= sigma_n; like -= tmp*tmp; sigma_n = like; sigma_s =scale;
			like = pow(like,-0.5f);
			yyy -= tmp; yyy *= like; // centering!
			mean = tmp;
			Tuple<double> proj;
			double bparam[3];
			double blike;
			// check for scale == 0;
			deriv[2] = 0.0f;
			yyy.show();
			param[0] = 0.0f; param[1] = 0.0f;
			double inner_product;
			/*
			do{
				scale = exp(param[1]);
				for(k=0,i=0;i<yyy.getSize();i++){
					for(j=0;j<i;j++) {covar.data[k++] = exp(param[0]);}
					covar.data[k++] = exp(param[0]) + scale;
				}
				like = (covar.Xformed_inner_product_of_inverse(yyy) + covar.log_determinant()) * -1.0f;
				proj = covar.divisionof(yyy);
				inner_product = proj[0] * proj[0];
				for(i=1;i<proj.getSize();i++) inner_product += proj[i] * proj[i];
				// trace ((KyyTK - K-1) * (K - scaleI))
				// trace ( KyyTK * K - I + scale * (K^-1  - KyyTK)  )
				//deriv[0] = (covar.trace_of_product(Trianglix<double>(proj)) - covar.getSize() + scale * (covar.trace_of_inverse() - inner_product) ) * exp(param[0]);
				deriv[0] = scale * covar.trace_of_inverse() - covar.getSize();
				printf("comps\t%e\t%e\t%e\t%e\t%e\n", covar.trace_of_product(Trianglix<double>(proj)) , (double) covar.getSize(), scale * covar.trace_of_inverse() , scale * inner_product,(double) covar.getSize() - scale * covar.trace_of_inverse());

				deriv[1] = inner_product - covar.trace_of_inverse();
				deriv[1] *= exp(param[1]);
				printf("ffcomps\t%e\t%e\n", covar.trace_of_inverse() ,  inner_product);
				// scale * trace ((KyyTK - K-1))

				memcpy(tparam,param,sizeof(double)*3);
				tparam[1]+= 0.0000001f;
				scale = exp(tparam[1]);
				for(k=0,i=0;i<yyy.getSize();i++){
					for(j=0;j<i;j++) {covar.data[k++] = exp(tparam[0]);}
					covar.data[k++] = exp(tparam[0]) + scale;
				}
				tlike = (covar.Xformed_inner_product_of_inverse(yyy) + covar.log_determinant()) * -1.0f;
				tparam[1]-= 0.0000002f;
				scale = exp(tparam[1]);
				for(k=0,i=0;i<yyy.getSize();i++){
					for(j=0;j<i;j++) {covar.data[k++] = exp(tparam[0]);}
					covar.data[k++] = exp(tparam[0]) + scale;
				}
				tlike -= (covar.Xformed_inner_product_of_inverse(yyy) + covar.log_determinant()) * -1.0f;


				printf("f[%e,%e] = %e   deriv = {%e,%e,%e}  {%e}\n", param[0], param[1], like,deriv[0],deriv[1],deriv[2], tlike * 5000000.0f);

				break;

			}while(!cs.updateAscent(like,param,deriv));*/

			blike = -(double)yyy.getSize(); // LL at 0,0,infinity
			bparam[0] = param[0];
			bparam[1] = param[1];
			bparam[2] = NAN;

			printf("start gg! f[0,0,inf] = %e\n", blike);

			for(unsigned int test=0;test<100;test++){
				cs.init(0.0000001f, 3);
				if (bparam[2] == NAN){
					param[0] = bparam[0] + 3.0f* sampleGaussian() - log(2.0f);
					param[1] = bparam[1] + 3.0f* sampleGaussian() - log(2.0f);
				}else{
					param[0] = bparam[0] + 0.5f* sampleGaussian() - log(2.0f);
					param[1] = bparam[1] + 0.5f* sampleGaussian() - log(2.0f);
				}
				param[2] = 10.0f * sampleGaussian() - log( 2.0f*sigma_s) + 2.0f* log((double)yyy.getSize());

				do{
					//scale = 1.0625f / (1.0f + exp(param[1]));
					scale = exp(param[1]);
					for(k=0,i=0;i<yyy.getSize();i++){
						for(j=0;j<i;j++) {tmp = vec[i].k - vec[j].k; covar.data[k++] = exp(param[0] - tmp*tmp * exp(param[2]));
							//printf("%e\t",exp(-tmp*tmp * exp(param[2])) );
						}
						covar.data[k++] = exp(param[0]) + scale;
					}//printf(" preLL=%e\n", like);
					like = (covar.Xformed_inner_product_of_inverse(yyy) + covar.log_determinant()) * -1.0f;
					if (!ExOp::isValid(like)||(param[2] < -2.0f -log( 2.0f*sigma_s))){if ((rand() & 255) != 0) test--; break;}
					proj = covar.divisionof(yyy);
					inner_product = proj[0] * proj[0];
					for(i=1;i<proj.getSize();i++) inner_product += proj[i] * proj[i];
					//deriv[0] = exp(param[0]) * covar.trace_of_product(Trianglix<double>(proj))  + scale * (covar.trace_of_inverse() - exp(param[0]) * inner_product) - (double)covar.getSize();
					deriv[0] = covar.trace_of_product(Trianglix<double>(proj))  + scale * (covar.trace_of_inverse() - inner_product) - (double)covar.getSize();

				//	printf("comps\t%e\t%e\t%e\t%e\t%e\n", covar.trace_of_product(Trianglix<double>(proj)) , (double) covar.getSize(), scale * covar.trace_of_inverse() , scale * inner_product,(double) covar.getSize() - scale * covar.trace_of_inverse());

					deriv[1] = inner_product - covar.trace_of_inverse();
					deriv[1] *= exp(param[1]);

					Trianglix<double> factrig = covar;
					for(k=0,i=0;i<yyy.getSize();i++){
						for(j=0;j<i;j++) {tmp = vec[i].k - vec[j].k; factrig.data[k++] *= -tmp*tmp;}
						factrig.data[k++] = 0.0f;
					}

					deriv[2] = factrig.trace_of_product(Trianglix<double>(proj)) - factrig.trace_of_division(covar);
					deriv[2] *= exp(param[2]);
					/*
					memcpy(tparam,param,sizeof(double)*3);
					tparam[0]+= 0.0000001f;
					scale = exp(tparam[1]);
					for(k=0,i=0;i<yyy.getSize();i++){
						for(j=0;j<i;j++) {tmp = vec[i].k - vec[j].k; covar.data[k++] = exp(tparam[0] - tmp*tmp * exp(tparam[2]));}
						covar.data[k++] = exp(tparam[0]) + scale;
					}
					tlike = (covar.Xformed_inner_product_of_inverse(yyy) + covar.log_determinant()) * -1.0f;
					tparam[0]-= 0.0000002f;
					scale = exp(tparam[1]);
					for(k=0,i=0;i<yyy.getSize();i++){
						for(j=0;j<i;j++) {tmp = vec[i].k - vec[j].k; covar.data[k++] = exp(tparam[0] - tmp*tmp * exp(tparam[2]));}
						covar.data[k++] = exp(tparam[0]) + scale;
					}
					tlike -= (covar.Xformed_inner_product_of_inverse(yyy) + covar.log_determinant()) * -1.0f;


					printf("f[%e,%e,%e] = %e   deriv = {%e,%e,%e}  %e\n", param[0], param[1], param[2], like,deriv[0],deriv[1],deriv[2], tlike * 5000000.0f);*/
				}while(!cs.updateAscent(like,param,deriv));

				if ((ExOp::isValid(like))&&(like > blike)){
					   printf("improved! f[%e,%e,%e] = %e (ite %i)\n", param[0], param[1], param[2], like, test);
					memcpy(bparam, param, sizeof(double)*3);
					blike = like;
				} else printf("rejected! f[%e,%e,%e] = %e (ite %i)\n", param[0], param[1], param[2], like, test);
			}
			sigma_s = exp(bparam[0]) * sigma_n;
			sigma_n *= 1.0625f / (1.0f + exp(bparam[1]));
			if (bparam[2] == NAN) {scale = 0.0f; sigma_s = 0.0f;}
			else scale = exp(-0.5f * param[2]);
		}
	}

// the inputs are SQUARED distances
double distanceToEllipse(double d1, double d2, double w, double f){
	double tmpb;
    //double distb;

	tmpb = d1 - d2;
	w *= 0.5f;
	w *= d1 + d2;
	f *= 0.5f;

	//distb = -2.0f * tmpb * tmpb;

	tmpb = d1 - d2;
	d1 += d2;
	w *= 0.5f;
	d1 = 4.0f * w * (d1+d2);
	d2 = -2.0f * tmpb * tmpb;


	f *= 0.5f;

/*	LFHPrimitive::PolyThing<double> poly;
	poly.setOrder(4);

	poly.coef[0] = w*w*(-8.0f*w*w+ 15*w*f-8.0f*f*f+d1) + (w-f)*(w-f)*d2;
	poly.coef[1] = 2.0f*(w*w*(-16.0f*w*w+ 23*w*f-8.0f*f*f+d1) + (w-f)*d2);
	poly.coef[2] = -48.0f*w*w+ 47*w*f-8.0f*f*f+ d1 + d2;
	poly.coef[3] = 16.0f * (f -2*w);
	poly.coef[4] = -8.0f;

	w =poly.newtonRoot(0.0f, f - w, LFHPrimitive::ExCo<double>::max());
	*/
	return(w);
}

double CubicRealRoot(double* coef, bool wantmiddle){
	double tmp = 3 *coef[3];
	double q = (coef[1] - (coef[2]*coef[2] / tmp) ) / tmp;
	double r = (-1.5f*coef[0] - ( ((coef[2]*coef[2]*coef[2] / tmp) - 1.5f*coef[2]*coef[1] )/tmp) )/tmp;
	tmp = q*q*q + r*r;
	//	printf("r:%f\td:%f\n",r,q);
	double s,t,ac,as;
	if ( tmp >= 0){ // only 1 real
		//		printf("preroot: s=%f\tt=%f\n",r + sqrt(q),r - sqrt(q));
		q = sqrt(tmp);
		tmp = r + q; s = pow(fabs(tmp),1/3.0f) * ((tmp <0) ? -1.0f:1.0f);
		tmp = r - q; t = pow(fabs(tmp),1/3.0f) * ((tmp <0) ? -1.0f:1.0f);
		//		printf("s=%f\tt=%f\n",s,t);
		return(s + t - coef[2] / (3*coef[3]));

    }else{ // all 3 are real, return the smallest one (closest to 0)
        tmp = r*r - tmp;
        t = acos(r/sqrt(tmp))/3;
        q = sqrt(-q); // q >= 0 is guarrantied!


        ac = cos(t);
        as = sin(t);
        if (wantmiddle){
            if (fabs(ac)*3.0f > fabs(as)) return((-ac +as* ((ac < 0.0f) ? sqrt(3.0f): -sqrt(3.0f) ))*q- coef[2] / (3*coef[3]));
            else return(ac*q*2 - coef[2] / (3*coef[3]));
        }else{

        s = ac*q*2 - coef[2] / (3*coef[3]);
        tmp = (-ac +as*sqrt(3.0f) )*q - coef[2] / (3*coef[3]);
        if (fabs(tmp) < fabs(s)) s =tmp;
        tmp = (-ac -as*sqrt(3.0f) )*q - coef[2] / (3*coef[3]);
        if (fabs(tmp) < fabs(s)) s =tmp;
        return(s);
        }
    }
}




double QuarticRealRoot(double* coef){
	double tmp = 3 *coef[4];
 //   double p =  (2.0f* coef[2] - (0.75f * coef[3]* coef[3] / coef[4])) / coef[4]; //(8*C2*C4-3*C3*C3)/(4*C4*C4)
	double q = (-4.0* coef[0] - ((3.0f * coef[3] * coef[1] - coef[2]*coef[2]) / tmp) ) / tmp;
	double r = ( 12.0f * coef[2] * coef[0] -4.5f * coef[1]*coef[1] + ((4.5f*coef[3]*(coef[2]*coef[1]-3.0f*coef[3]*coef[0]) -coef[2]*coef[2]*coef[2] )/tmp) )/ (tmp*tmp);
	// =(24*C2*C0-9*(C1*C1)+(9*C3*(C2*C1-3*C3*C0)-2*C2*C2*C2)/(3*C4) ))/(18*C4*C4)
    tmp = q*q*q + r*r;
    double root,t,s,ac,as;
    if (tmp > 0) { // one real root for cubic
        // since imaginary roots are roots have equal opposite imaginary components, so they describe a real polynomial, the real root correspond to the
        // factorization that pairs properly imaganinary roots.
		q = sqrt(tmp);
		tmp = r + q; s = pow(fabs(tmp),1/3.0f) * ((tmp <0) ? -1.0f:1.0f);
		tmp = r - q; t = pow(fabs(tmp),1/3.0f) * ((tmp <0) ? -1.0f:1.0f);
		//		printf("s=%f\tt=%f\n",s,t);
		root = s + t - coef[2] / (3*coef[3]);
    }else{ // three real roots for cubic
		tmp = r*r - tmp;
		t = acos(r/sqrt(tmp))/3;
		q = sqrt(-q); // q >= 0 is guarrantied!

		ac = cos(t);
		as = sin(t);

		s = ac*q*2 - coef[2] / (3*coef[3]);
		tmp = (-ac +as*sqrt(3.0f) )*q - coef[2] / (3*coef[3]);
		if (fabs(tmp) < fabs(s)) s =tmp;
		tmp = (-ac -as*sqrt(3.0f) )*q - coef[2] / (3*coef[3]);
		if (fabs(tmp) < fabs(s)) s =tmp;
		root =  s;
    }
    printf("cubic root: %f\n",root);

    // discriminent!
    tmp = -0.25f * root + ( (coef[2] + (-0.375f * coef[3] * coef[3] / coef[4])) - (coef[1] + 0.5f*coef[3] * ((0.25f * coef[3]* coef[3] / coef[4]) - coef[2]) / coef[4]) * pow(fabs(root), -0.5f) ) / ( 2.0f * coef[4]);
    if (tmp < 0.0f) return(-1.0f);
    // should *not* be negative!
    r = 0.25f * coef[3] / coef[4] + sqrt(fabs(root));
    return (sqrt(tmp) > r) ? r + sqrt(tmp) : r - sqrt(tmp);
}



// point at -3,-1,1,3 return closest to center, valid root
double CubicInterpolationRoot(double* pt , bool is_monotonic){

	double coef[4];
	double q,r;
	q = (pt[3] - pt[0])/ (3 * 16.0f);
	r = (pt[2] - pt[1])/16.0f;
	coef[1] = -q + 9.0f *r;
	coef[3] = q - r;
	q = (pt[3] + pt[0])/ 16.0f;
	r = (pt[2] + pt[1]) /(16.0f);
	coef[0] = -q + 9.0f *r;
	coef[2] = q - r;



	double tmp = 3 *coef[3];
	q = (coef[1] - (coef[2]*coef[2] / tmp) ) / tmp;
	r = (-1.5f*coef[0] - ( ((coef[2]*coef[2]*coef[2] / tmp) - 1.5f*coef[2]*coef[1] )/tmp) )/tmp;
	tmp = q*q*q + r*r;
	//	printf("r:%f\td:%f\n",r,q);
	double s,t,ac,as;
	if ( tmp >= 0){ // only 1 real
		//		printf("preroot: s=%f\tt=%f\n",r + sqrt(q),r - sqrt(q));
		q = sqrt(tmp);
		tmp = r + q; s = pow(fabs(tmp),1/3.0f) * ((tmp <0) ? -1.0f:1.0f);
		tmp = r - q; t = pow(fabs(tmp),1/3.0f) * ((tmp <0) ? -1.0f:1.0f);
		//		printf("s=%f\tt=%f\n",s,t);
		return(s + t - coef[2] / (3*coef[3]));

	}else{ // all 3 are real, return the smallest one (closest to 0)
		tmp = r*r - tmp;
		t = acos(r/sqrt(tmp))/3;
		q = sqrt(-q); // q >= 0 is guarrantied!


		ac = cos(t);
		as = sin(t);

//		if (is_monotonic){
			// ignores illegal roots!
//		}else{
			s = ac*q*2 - coef[2] / (3*coef[3]);
			tmp = (-ac +as*sqrt(3.0f) )*q - coef[2] / (3*coef[3]);
			if (fabs(tmp) < fabs(s)) s =tmp;
			tmp = (-ac -as*sqrt(3.0f) )*q - coef[2] / (3*coef[3]);
			if (fabs(tmp) < fabs(s)) s =tmp;
		return(s);
//		}



	}
}

double sinc(double xx){

	if (fabs(xx) < pow(0.5f,23.0f)) return(1.0f);
	else return(sin(xx) / xx);

	}

double sinc_d(double xx){
	// (sin x - x cos x) / x^2

		if (fabs(xx) < pow(0.5f,23.0f)) return((-2.0f  / 3.0f) * xx);
	else return(((sin(xx) / xx) - cos(xx))/ xx);

	}




SetComparison::SetComparison(){}
SetComparison::SetComparison(int v) : comp(v){}

WarH static_warning_handdle;

void ArgumentParser::readList(char* const text, vector<unsigned int> & _out){

	char* cur = text;
	bool range = false;
	unsigned int tmpint;
	unsigned int tmpint_ite;
	while(true){
		switch(*cur){
			case '\0': return;
			case ' ': break;
			case '-': range = true;
			default:
				if (((*cur) >= '0')&&((*cur) <= '9')){
					tmpint = atoi(cur);
					while(((*cur) >= '0')&&((*cur) <= '9')) cur++;
					cur--;
					if (range){
						tmpint_ite =  _out[_out.size()-1];
						if (tmpint >= tmpint_ite) for(;tmpint_ite <=tmpint ;tmpint_ite++) _out.push_back(tmpint_ite);
						else for(;tmpint_ite >=tmpint ;tmpint_ite--) _out.push_back(tmpint_ite);
					}else _out.push_back(tmpint);
			}	}
			cur++;
	}	}
	void ArgumentParser::readList(char* const text, Vector<unsigned int> & _out){

		char* cur = text;
		bool range = false;
		unsigned int tmpint;
		unsigned int tmpint_ite;
		while(true){
			switch(*cur){
				case '\0': return;
				case ' ': break;
				case '-': range = true;
				default:
				if (((*cur) >= '0')&&((*cur) <= '9')){
					tmpint = atoi(cur);
					while(((*cur) >= '0')&&((*cur) <= '9')) cur++;
					cur--;
					if (range){
						tmpint_ite =  _out[_out.getSize()-1];
						if (tmpint >= tmpint_ite) for(;tmpint_ite <=tmpint ;tmpint_ite++) _out.push_back(tmpint_ite);
						else for(;tmpint_ite >=tmpint ;tmpint_ite--) _out.push_back(tmpint_ite);
					}else _out.push_back(tmpint);
			}	}
			cur++;
	}	}

void ArgumentParser::help_routine(char c){
	switch(c){
		case 'l': printf("* Copyright (C) 2013 Louis-Francois Handfield\n* e-mail: lfhandfield@gmail.com\n* \n* This program is free software; uppon the notification to the licensor\n* of the licencee identity and nature of use, the licencee can redistribute this\n* program and/or modify it under the terms of the GNU General Public\n* License as published by the Free Software Foundation; either version 2\n* of the License, or (at the licencee option) any later version. As such,\n* no further notification of any kind are required.\n*\n* This program is distributed in the hope that it will be useful, but\n* WITHOUT ANY WARRANTY; without even the implied warranty of\n* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU\n* General Public License for more details.\n*\n* You should have received a copy of the GNU General Public License\n* along with this program; if not, write to the Free Software\n* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.\n");
		case 's': printf("Optionnal Argument syntax:\n\nFlags are separated by space only. Tab character are ignored.\n\nExample 1:\n\t-a[BC] (file x) (file y)\nThe flag 'a' requires 2 arguments, which are 2 file path. B and C are flag specific options, which are letters appended to the flag (no space, any order). The option are described below the flag. Here is some valid used of the flag:\n\t(program) -aB circle.tif output.txt\n\t(program) -aBC squarre.tif output.txt\n\nExample 2:\n\t-b (int x) (int y = 0)\n The flag b has no options, and takes 1 or 2 arguments. If a single argument is provided, the parameter will match the value shown after the '=' sign. Example:\n\t(program) -b 2 4\n\t(program) -b 7\t\t(which is equivalent to '(program) -b 7 0' )\n\nSpecial flag:\n a lone '-' is a special flag, which has 2 possible meanings:\n First, it can be used to before default argument to the process; allows to place such argument anywhere, and also allows to resolve term when both default and the last flag uses variable number of arguments.\n Second, if it is used as the very first flag, it takes one arguement, which is a file path to a file containing flags for the process; flags within the file will be read, so it would be equivalent to substitute this custom flag by the file content.\n"); break;
		default:
			help();
			printf("\t-hL: for source code licence for this software\n");
			printf("\t-hS: for help of flag syntax\n");
			printf("Compiled on the %s, at %s\n",__DATE__, __TIME__);
	}
	LFH_ALIVE; exit(1);
	}

int ArgumentParser::operator()(const char * argslist){
    Vector<char*> parsedargs;
    const char* st = argslist;
    unsigned int length=0;
    bool backsl = false;
    bool inquote = false;
    parsedargs.push_back(NULL);
    unsigned int i;
    for(const char* cur = argslist;*cur != '\0';cur++){
        switch(*cur){
            case '\\':
                    if (inquote) length++;
                    else{
                    if (backsl) length++;
                    backsl = !backsl;
                    }
                break;
            case '\'':
                    if (backsl) {length++; backsl =false;}
                    else inquote = !inquote;
                break;
            case ' ':
                if (inquote) length++;
                else if (backsl) {length++; backsl = false;}
                else if (length != 0){
                    parsedargs.push_back(new char[length+1]);
                    parsedargs.last()[length] = '\0';
                    for(cur=st, i=0;i < length;cur++){
                        switch(*cur){
                            case '\\':
                            if (inquote) parsedargs.last()[i++] = *cur;
                            else{
                            if (backsl) parsedargs.last()[i++] = *cur;
                            backsl = !backsl;
                            }
                            break;
                            case '\'':
                                if (backsl) {parsedargs.last()[i++] = *cur; backsl =false;}
                                else inquote = !inquote;
                            break;
                            case ' ':
                                parsedargs.last()[i++] = *cur;
                                backsl = !backsl;
                            break;
                            default:
                               parsedargs.last()[i++] = *cur;
                        }
                    }
                    length = 0;
                    if (inquote) {cur++; inquote = false;}
                    while(*cur == ' ') cur++;
                    st = cur--;backsl =false;
                }else st = cur+1;

            break;
            default:
            length++;
            }
        }
    if (length != 0){
        parsedargs.push_back(new char[length+1]);
        parsedargs.last()[length] = '\0';i=0;
        for(const char* cur=st;i < length;cur++){
            switch(*cur){
                case '\\':
                    if (inquote) parsedargs.last()[i++] = *cur;
                    else{
                    if (backsl) parsedargs.last()[i++] = *cur;
                    backsl = !backsl;
                    }
                break;
                case '\'':
                    if (backsl) {parsedargs.last()[i++] = *cur; backsl =false;}
                    else inquote = !inquote;
                break;
                case ' ':
                    parsedargs.last()[i++] = *cur;
                    backsl = !backsl;
                break;
                default:
                   parsedargs.last()[i++] = *cur;
            }
        }
    }

    return (*this)(parsedargs.getSize(), &(parsedargs[0]));
}

int ArgumentParser::operator()(int nbargs, char * const * args, bool has_prog_name_as_argument){
	int min, max,i,j;
	if (has_prog_name_as_argument) nbargs--;
	char * const * cargs = args +((has_prog_name_as_argument) ? 1 : 0);
	int min_default;

	char emptyfun[] = "";
	max =0;nbaddtoken(emptyfun, min_default, max);


	while((nbargs > min_default)&&(cargs[0][0] == '-')){

		if (cargs[0][1] == '\0') {nbargs--; cargs++; break;}
		else if (cargs[0][1] == 'h') help_routine(cargs[0][2]);

		max =0;nbaddtoken((*cargs) + 1, min, max);

		if ((min == max)&&(max != 0)){
			store(cargs,nbargs);
			return defstore(cargs,nbargs);
		}else{
			if (max < min) max =min;
			for(i=1;i<nbargs;i++) {
				if (i >= nbargs - min_default) break;
				if (cargs[i][0] == '-') { //number? ignore
					for(j=1; cargs[i][j] == ' ';j++);
					if (((cargs[i][j] < '0')||(cargs[i][j] > '9'))&&(cargs[i][j] != '.')) break;
				}
			}
			if ( i-1 > max){
				if (i >= nbargs - min_default) i = max+1;
				else{
				if (max > min) fprintf(stderr,"flag '%s' needs %i to %i arguemnts, found %i\n",(*cargs) + 1 , min, max,i-1);
				else fprintf(stderr,"flag '%s' needs %i arguments, found %i\n", (*cargs) + 1 , min,i-1);
				fprintf(stderr,"for more details, try %s -h\n", *args);
				exit(1);
				}
			}
			if (i-1 < min) {
				if (max > min) fprintf(stderr,"flag '%s' needs %i to %i arguemnts, found %i\n",(*cargs) + 1 , min, max,i-1);
				else fprintf(stderr,"flag '%s' needs %i arguments, found %i\n", (*cargs) + 1 , min,i-1);
				fprintf(stderr,"for more details, try %s -h\n", *args);
				exit(1);
			}
			store(cargs,i-1);
			cargs += i;
			nbargs -= i;
	}	}

	for(i=0;i<nbargs;i++){
		if (cargs[i][0] == '-'){
			if (cargs[i][1] == '\0') {nbargs--; cargs++; break;}
			else if (cargs[i][1] == 'h') help_routine(cargs[0][2]);
			max =0;	nbaddtoken((*cargs) + 1, min, max);
			if ((min == max)&&(max != 0)){
				store(cargs,nbargs);
				return defstore(cargs,nbargs);
	}	}	}
	// no flag = default flag \0, (last flag too)
	max =0; nbaddtoken(emptyfun, min, max);
	if (max < min) max =min;
	if ((nbargs < min)||(nbargs > max)){
		if (max > min) fprintf(stderr,"process needs %i to %i arguemnts, found %i\n", min, max,nbargs);
		else fprintf(stderr,"process needs %i arguments, found %i\n", min,nbargs);
		fprintf(stderr,"for more details, try %s -h\n", *args);
		exit(1);
	}
	return defstore(cargs,nbargs);
}

	TableReader::TableReader(const char* path, const Vector<const char *> &_columns){
		Vector< KeyElem<unsigned int, unsigned int> > loc_ibo;
		unsigned int i,j,k,l;
		myHashmap<string, unsigned int> all_needed_colname;
		char buffer[65536];
		char sep;
		float tmpfl;
		unsigned int extra_read =0; unsigned int found=0;
		for(j=0;j<_columns.size();j++) {
			if ((_columns[j][0] >= '0')&&(_columns[j][0] <= '9')){ // check if whole string in integer, so it is interpreted as a column index
				for(i= strlen(_columns[j])-1;i!=0;i--) if ((_columns[j][i] < '0')||(_columns[j][i] > '9')) break;
			} else i = 1;
			if (i ==0) loc_ibo.push_back(KeyElem<unsigned int, unsigned int>(atoi(_columns[j]), j));
			else{
				k=0;
				for(i=0;_columns[j][i] != '\0';i++){
					if ((_columns[j][i] == '-') || (_columns[j][i] == '+')|| (_columns[j][i] == '/')|| (_columns[j][i] == '*')|| (_columns[j][i] == '(')|| (_columns[j][i] == ')')|| (_columns[j][i] == ',')){
						memcpy(buffer, _columns[j] + k, i-k);
						if ((buffer[0] >= '0')&&(buffer[0] <= '9')){ // number!
							tmpfl = atof(buffer); constants.push_back(tmpfl);
							// if (_columns[j][i] == '(') implicit multiplication... todo
						}else{
							buffer[i-k] = '\0';
							if (_columns[j][i] == '('){
								if ((strcmp(buffer,"log") == 0)||(strcmp(buffer,"log10") == 0)||(strcmp(buffer,"log2") == 0)||(strcmp(buffer,"exp") == 0)||(strcmp(buffer,"exp10") == 0)||(strcmp(buffer,"exp2") == 0)||(strcmp(buffer,"sin") == 0)||(strcmp(buffer,"cos") == 0)||(strcmp(buffer,"sqrt") == 0)||(strcmp(buffer,"pow") == 0)||(strcmp(buffer,"pi") == 0)){
									k=i+1;
								} else {
									l = all_needed_colname.find(string(buffer));
									if (l == 0xFFFFFFFF) {all_needed_colname[string(buffer)] = 0xFFFFFFFF; extra_read++;}
									k=i+1;
								}
							}else if (i != k){
							l = all_needed_colname.find(string(buffer));
							if (l == 0xFFFFFFFF) {all_needed_colname[string(buffer)] = 0xFFFFFFFF; extra_read++;}

							k=i+1;
							}else k++;
						}
					}
				}

				if ((_columns[j][k] >= '0')&&(_columns[j][k] <= '9')){
					tmpfl = atof(buffer); constants.push_back(tmpfl);
				}else{

					if (k == 0){
						l = all_needed_colname.find(string(_columns[j]));
						if (l == 0xFFFFFFFF) all_needed_colname[string(_columns[j])] = j;
						else if (all_needed_colname.deref(l) == 0xFFFFFFFF) {all_needed_colname.deref(l) = j; extra_read--;}
						else found++; // needs copying

					}else{
						if (_columns[j][k] != '\0'){
						l = all_needed_colname.find(string(_columns[j] + k));
						if (l == 0xFFFFFFFF) {all_needed_colname[string(_columns[j] + k)] = 0xFFFFFFFF;extra_read++;}
						}
						found++;
					}
				}
			}
		}
	//	printf("gothere!\n"); fflush(stdout);
//ExOp::show(all_needed_colname);
		loc_ibo.sort();


		f = fopen(path, "r+"); if (f ==NULL) {fprintf(stderr, "Could not open file %s!\n", path); exit(1);}
		nbcols = _columns.size(); nbread =nbcols + extra_read - found;
	//	printf("special columns: %i, extra reads: %i nb columns: %i\n",found,extra_read,nbcols );
		read_info = new unsigned int[nbread*2];
		col_name = new string[nbcols + extra_read];
	//	printf("%i\t%i\t%i\n", nbcols, extra_read, found); fflush(stdout);
		// fill col_name
		unsigned int cur_num_index=0;

		for(i=0;i<_columns.size();i++){
			if ((cur_num_index >= loc_ibo.size())||(loc_ibo[cur_num_index].k != k)) col_name[i] = string(_columns[i]);
			else col_name[i] = string("+"); // imposible to match
		}
		j = _columns.size();
		for(i=0;i<all_needed_colname.getSize();i++){
			if (all_needed_colname.deref(i) == 0xFFFFFFFF){
				col_name[j] = all_needed_colname.deref_key(i);
				all_needed_colname.deref(i) = j++;
			}
		}

		found=0;

		read_info[0] =0;k=0;cur_num_index=0;j=0;
		while(2 == fscanf(f,"%[^\t\n\r]%c",buffer, &sep)){
			if ((cur_num_index < loc_ibo.size())&&(loc_ibo[cur_num_index].k == k)) {i = loc_ibo[cur_num_index++].d;col_name[found] = string(buffer);}
			else for(i=0;i<nbcols + extra_read;i++) if (strcmp(buffer, col_name[i].c_str()) == 0) break;

			if (i == nbcols + extra_read) read_info[(found << 1)]++;
			else{
				read_info[(found << 1) | 1] = i;
				found++;

				if (found == nbread) break;
				read_info[(found << 1)] =0;
			}
			k++;

			if ((sep == '\n')||(sep == '\r')) break;
		}
		if (found != nbread){
			// columns are missing, report em!
			for(i=0;i<found;i++) col_name[read_info[(i <<1) | 1]].clear();
			fprintf(stderr,"tableread error: missing %i columns out of %i! The missing column name are ammong the following:\n",  nbread - found , nbread );
			for(i=0;i<nbcols + extra_read;i++) if (col_name[i].length() != 0) fprintf(stderr,"'%s'\n", col_name[i].c_str());
			exit(1);
		}
		if (j < loc_ibo.size()) {fprintf(stderr,"table file %s does not a %ith column!\n", path, loc_ibo[j].k);exit(1);}

		while((sep != '\n')&&(sep != '\r'))	fscanf(f,"%[^\t\n\r]%c",buffer, &sep); // flush the resto of the row

		current_row = new string[nbcols + extra_read];

		// operattions!

		Vector<unsigned int> delay;
		bool isterm;

		for(j=0;j<_columns.size();j++) {
			if ((_columns[j][0] >= '0')&&(_columns[j][0] <= '9')){
				for(i= strlen(_columns[j])-1;i!=0;i--) if ((_columns[j][i] < '0')||(_columns[j][i] > '9')) break;
			} else i = 1;
			if (i !=0){
				k=0;sep = '\0';
				for(i=0;_columns[j][i] != '\0';i++){
					if ((_columns[j][i] == '-') || (_columns[j][i] == '+')|| (_columns[j][i] == '/')|| (_columns[j][i] == '*')|| (_columns[j][i] == '(')|| (_columns[j][i] == ')')|| (_columns[j][i] == ',')){
						memcpy(buffer, _columns[j] + k, i-k);
						buffer[i-k] = '\0';
					//	printf("%s\n", buffer);fflush(stdout);
						isterm = false;
						if ((_columns[j][k] >= '0')&&(_columns[j][k] <= '9')){
							oper.push_back(0x80000010);
							if (_columns[j][i] == '(') {sep = '*';}
							else isterm = true;
						}else if (_columns[j][i] == '('){
							switch(buffer[0]){
								case 0:delay.push_back(0x80000000);break; // simple parenthesis
								case 'e':
									if (strcmp(buffer, "exp") == 0) delay.push_back(0x80000006);
									else if (strcmp(buffer, "exp2") == 0) delay.push_back(0x80000008);
									else if (strcmp(buffer, "exp10") == 0) delay.push_back(0x8000000A);
									else {oper.push_back(all_needed_colname[string(buffer)]); sep = '*';}
									break;
								case 'l':
									if (strcmp(buffer, "log") == 0) delay.push_back(0x80000005);
									else if (strcmp(buffer, "log2") == 0) delay.push_back(0x80000007);
									else if (strcmp(buffer, "log10") == 0) delay.push_back(0x80000009);
									else {oper.push_back(all_needed_colname[string(buffer)]); sep = '*';}
								break;
								case 's':
									if (strcmp(buffer, "sin") == 0) delay.push_back(0x8000000B);
									else if (strcmp(buffer, "sqrt") == 0) delay.push_back(0x8000000E);
									else {oper.push_back(all_needed_colname[string(buffer)]); sep = '*';}
								break;
								case 'p':
									if (strcmp(buffer, "pow") == 0) delay.push_back(0x8000000F);
									else if (strcmp(buffer, "pi") == 0) {oper.push_back(0x80000011); isterm = true; if (_columns[j][++i] != ')') {fprintf(stderr, "tableread error: pi() takes no argument\n"); exit(1);}i++;}
									break;
								default:
									if (strcmp(buffer, "cos") == 0) delay.push_back(0x8000000C);
									else if (strcmp(buffer, "rad") == 0) delay.push_back(0x8000000D);
									else {oper.push_back(all_needed_colname[string(buffer)]); sep = '*';}
							}

						}else {if  (i != k) oper.push_back(all_needed_colname[string(buffer)]); isterm = true;}

						if (isterm){
						switch(sep){
							case '+': delay.push_back(0x80000001);break;
							case '-': delay.push_back(0x80000002);break;
							case '*': oper.push_back(0x80000003);break;
							case '/': oper.push_back(0x80000004);break;
							case ',': delay.push_back(0x80000012);break;
						}

						if (_columns[j][i] == ')'){
			//				printf("HALLEO1\n");
			//				printf("%i\n",delay.size()); fflush(stdout);
			//				printf("%i\n",delay[delay.size()-1]);fflush(stdout);
							do{
							if ((delay[delay.size()-1] != 0x80000000)&&(delay[delay.size()-1] != 0x80000012)) oper.push_back(delay[delay.size()-1]);
							}while((oper[oper.size()-1] == 0x80000001)||(oper[oper.size()-1] == 0x80000002));
							delay.pop_back();
							sep = '\0';
						}else sep = _columns[j][i];
						}
						k = i+1;
					}
				}

//

				if (k != 0){
	//				printf("operations!:\n");for(i=0;i<oper.size();i++) printf("%X\n", oper[i]);
	//				printf("delayed!:\n");for(i=0;i<delay.size();i++) printf("%X\n", delay[i]);
					if (_columns[j][k] != ')'){
					if (_columns[j][k] != '\0') oper.push_back(all_needed_colname[string(_columns[j] + k)]);
					switch(sep){
						case '+': oper.push_back(0x80000001);break;
						case '-': oper.push_back(0x80000002);break;
						case '*': oper.push_back(0x80000003);break;
						case '/': oper.push_back(0x80000004);break;
					}

					}
					for(i=delay.size();i >0 ;) {if (delay[--i] != 0x80000000) oper.push_back(delay[i]);}
					delay.toMemfree();


					oper.push_back(0x40000000 | j);
				}else{
					if (all_needed_colname[string(_columns[j])] != j){ // multiple uses... copy!
						oper.push_back(all_needed_colname[string(_columns[j])]);
						oper.push_back(0x40000000 | j);
					}
				}
			}
		}


	//	printf("operations!:\n");for(i=0;i<oper.size();i++) printf("%X\n", oper[i]);
	//	fflush(stdout);
	}
	TableReader::~TableReader(){fclose(f);
		delete[](read_info);
		delete[](current_row);
		delete[](col_name);
	}

	bool TableReader::nextRow(){
		double dbuffer[16];
		unsigned int i,j;
		char buffer[65536];
		char sep;
		unsigned int constant_offset;
		string current;
		for(j=0;j<nbread;j++){
			for(i=0;i<=read_info[(j << 1)];i++) {
				if (2 != fscanf(f,"%[^\t\n\r]%c",buffer, &sep)){
					if (1 == fscanf(f,"%c", &sep)){
						buffer[0] = '\0';
						if ((sep == '\n')||(sep == '\r')) break;
					}else return false;
				}
			}
			if (i<read_info[(j << 1)]) return false;
			current_row[read_info[(j << 1) | 1]] = string(buffer);
			if ((i == read_info[(j << 1)])&&(j != nbcols-1)) return false;
		}
		while((sep != '\n')&&(sep != '\r'))	if (2 != fscanf(f,"%[^\t\n\r]%c",buffer, &sep)) fscanf(f,"%c", &sep);


		i=0;constant_offset=0;
		for(j=0;j<oper.size();j++){

			if (oper[j] &0x80000000){
				switch(oper[j] &0x7FFFFFFF){
					case 1: i--;dbuffer[i-1] += dbuffer[i]; break;
					case 2: i--;dbuffer[i-1] -= dbuffer[i]; break;
					case 3: i--;dbuffer[i-1] *= dbuffer[i]; break;
					case 4: i--;dbuffer[i-1] /= dbuffer[i]; break;
					case 5: dbuffer[i-1] = log(dbuffer[i-1]); break;
					case 6: dbuffer[i-1] = exp(dbuffer[i-1]); break;
					case 7: dbuffer[i-1] = log2(dbuffer[i-1]); break;
					case 8: dbuffer[i-1] = exp2(dbuffer[i-1]); break;
					case 9: dbuffer[i-1] = log10(dbuffer[i-1]); break;
					case 10: dbuffer[i-1] = exp(dbuffer[i-1] * log(10.0f)); break;
					case 11: dbuffer[i-1] = sin(dbuffer[i-1]); break;
					case 12: dbuffer[i-1] = cos(dbuffer[i-1]); break;
					case 13: dbuffer[i-1] *= 2.0f * M_PI; break;
					case 14: dbuffer[i-1] = sqrt(dbuffer[i-1]); break;
					case 15: i--;dbuffer[i-1] = pow(dbuffer[i-1], dbuffer[i]); break;
					case 16: dbuffer[i++] = constants[constant_offset++]; break;
					case 17: dbuffer[i++] = M_PI; break;
					case 18: buffer[0] = ','; sprintf(buffer+1, "%.32e", dbuffer[i]); current = string(buffer) + current; break; // concatenate!
				}
			}else if (oper[j] &0x40000000){ // write
				i--;//printf("w\n");fflush(stdout);
				sprintf(buffer, "%.32e", dbuffer[i]);
				//printf("%s!\n",buffer);
				current_row[(oper[j] &0x3FFFFFFF)] = string(buffer) + current;
				current.clear();
			}else{//printf("r\n");fflush(stdout);
				dbuffer[i++] = atof(current_row[oper[j]].c_str());
			}
		}
	//	printf("end\n");fflush(stdout);
		return true;
	}

	TiffFile::TiffFile(const char* path, bool writeonly) : curfp_pos(4), endfile_pos(0){
		if (path != NULL){
			f = (writeonly) ? NULL : fopen(path,"rb+");
	char buffer[256];

	if (f == NULL){ // it's a new file, init header!
		f = fopen(path,"wb+");
		if (f ==NULL) {fprintf(stderr,"could not open/create %s!\n", path); return;}
		buffer[0] = 'I';
		buffer[1] = 'I';
		*(short*)(buffer + 2) = 42;
		*(short*)(buffer + 2) = 42;
		*(int*)(buffer + 4) = 0;
		fwrite(buffer,sizeof(char),8,f);

	}else{ // filed exists!
		fread(buffer,sizeof(char),4,f);
		if (buffer[0] == 'M'){
			if (buffer[1] != 'M') {fprintf(stderr,"tried to open %s, with is not a tiff file!\n", path); exit(1);}
			inv =true;
			if (*((short*)(buffer + 2)) != 10752) {fprintf(stderr,"tried to open %s, with is not a tiff file!\n", path); exit(1);}
			} else if (buffer[0]  == 'I'){
			if (buffer[1] != 'I') {fprintf(stderr,"tried to open %s, with is not a tiff file!\n", path); exit(1);}
			inv=false;
		if (*((short*)(buffer + 2)) != 42) {fprintf(stderr,"tried to open %s, with is not a tiff file!\n", path); exit(1);}
		}else fprintf(stderr,"tried to open %s, with is not a tiff file!\n", path);

	}
	} else f = NULL;
}

TiffFile::~TiffFile(){
	if (f != NULL) fclose(f);
}

bool TiffFile::gotoNext(){
	unsigned int i;
	unsigned short s;
	fseek(f, curfp_pos, SEEK_SET);
	fread(&i,sizeof(unsigned int),1,f);
	if (inv) ExOp::bytereverse(i);
	if (i == 0) {fseek(f, -4 * sizeof(char), SEEK_CUR); return(false);}
	fseek(f, i, SEEK_SET);
	fread(&s,sizeof(unsigned short),1,f);
	if (inv) ExOp::bytereverse(s);
	curflaglist.setSize(s*12);
	fread(&(curflaglist[0]),sizeof(unsigned char),s*12,f);
	curfp_pos = i + (2+s*12)*sizeof(unsigned char);
	return(true);
}

int TiffFile::flagType(int flagindex){
		unsigned short type = (*(unsigned short*)&(curflaglist[flagindex*12+2]));
		if (inv) ExOp::bytereverse(type);
		int out;
		switch(type){
			case 1:
			case 2:
			case 6:
			case 7:
				out = (*(unsigned int*)&(curflaglist[flagindex*12+4]));
				if (inv) ExOp::bytereverse(out);
				if (out <= 4) out =  type | 0x00000100;
				else out = type;
			break;
			case 3:
			case 8:
				out = (*(unsigned int*)&(curflaglist[flagindex*12+4]));
				if (inv) ExOp::bytereverse(out);
				if (out <= 2) out =  type | 0x00000100;
				else out = type;
			break;
			case 4: // int
			case 9:
			case 11:
				out = (*(unsigned int*)&(curflaglist[flagindex*12+4]));
				if (inv) ExOp::bytereverse(out);
				if (out == 1) out =  type | 0x00000100;
				else out = type;
			break;
		}
	return(out);
	}

template< >
int TiffFile::getValue<int>(int flagindex){
	int type = flagType(flagindex);
	switch(type){
		case 0x00000004:
			fseek(f, *(unsigned int*)&(curflaglist[flagindex*12+ 8]), SEEK_SET);
		break;
		case 0x00000101:
		case 0x00000102: {unsigned char c = *(unsigned char*)&(curflaglist[flagindex*12+ 8]); if (inv) ExOp::bytereverse(c); return(c);}
		case 0x00000103: {unsigned short c = *(unsigned short*)&(curflaglist[flagindex*12+ 8]); if (inv) ExOp::bytereverse(c); return(c);}
		case 0x00000104: {unsigned int c = *(unsigned int*)&(curflaglist[flagindex*12+ 8]); if (inv) ExOp::bytereverse(c); return(c);}
		case 0x00000106: {char c = *(char*)&(curflaglist[flagindex*12+8]); if (inv) ExOp::bytereverse(c); return(c);}
		case 0x00000108: {short c = *(short*)&(curflaglist[flagindex*12+8]); if (inv) ExOp::bytereverse(c); return(c);}
		case 0x00000109: {int c = *(int*)&(curflaglist[flagindex*12+ 8]); if (inv) ExOp::bytereverse(c); return(c);}
	}
	return(0);
}

template< >
unsigned int TiffFile::getValue<unsigned int>(int flagindex){
		int type = flagType(flagindex);
		int rout;
		unsigned int i_out;
		switch(type){
			case 0x00000003: {//:
				unsigned short i_buf;
				fseek(f, *(unsigned int*)&(curflaglist[flagindex*12+8]), SEEK_SET);
				if ((rout = fread(&i_buf,sizeof(unsigned short), 1, f)) != 1) return 0;
				i_out =i_buf;
				}//:
				return(i_out);
			case 0x00000004:
				fseek(f, *(unsigned int*)&(curflaglist[flagindex*12+8]), SEEK_SET);
				if ((rout = fread(&i_out,sizeof(unsigned int), 1, f)) != 1) return 0;
				return(i_out);
			case 0x00000101:
		case 0x00000102: {unsigned char c = *(unsigned char*)&(curflaglist[flagindex*12+8]); if (inv) ExOp::bytereverse(c); return(c);}
		case 0x00000103: {unsigned short c = *(unsigned short*)&(curflaglist[flagindex*12+8]); if (inv) ExOp::bytereverse(c); return(c);}
		case 0x00000104: {unsigned int c = *(unsigned int*)&(curflaglist[flagindex*12+8]); if (inv) ExOp::bytereverse(c); return(c);}
		case 0x00000106: {char c = *(char*)&(curflaglist[flagindex*12+8]); if (inv) ExOp::bytereverse(c); return(c);}
		case 0x00000108: {short c = *(short*)&(curflaglist[flagindex*12+8]); if (inv) ExOp::bytereverse(c); return(c);}
		case 0x00000109: {int c = *(int*)&(curflaglist[flagindex*12+8]); if (inv) ExOp::bytereverse(c); return(c);}
		}
		return(0);
	}

	template< >
	vector<unsigned int> TiffFile::getValue< vector< unsigned int> >(int flagindex){
		int type = flagType(flagindex);
        int rout;
		vector< unsigned int> f_out;
		unsigned int tmptmp;
		unsigned short tmptmps;
		unsigned int i = *(unsigned int*)&(curflaglist[flagindex*12+8]);
		unsigned int arraysize = (*(unsigned int*)&(curflaglist[flagindex*12+4]));
		if (inv) {ExOp::bytereverse(i);ExOp::bytereverse(arraysize);}
		switch(type){
			case 0x00000004:
				fseek(f, i, SEEK_SET);
				for(i=0;i< arraysize ;i++) {
					if ((rout = fread(&tmptmp,sizeof(unsigned int), 1, f)) != 1) return f_out;
					if (inv) ExOp::bytereverse(tmptmp);
					f_out.push_back(tmptmp);
				}
				break;
			case 0x00000003:
				fseek(f, i, SEEK_SET);
				for(i=0;i< arraysize ;i++) {
					if ((rout = fread(&tmptmps,sizeof(unsigned short), 1, f)) != 1) return f_out;
					if (inv) ExOp::bytereverse(tmptmps);
					f_out.push_back(tmptmps);
				}
				break;
			case 0x00000104:
				f_out.push_back(i);
				break;
			case 0x00000103:
				tmptmps = *(unsigned short*)&(curflaglist[flagindex*12+8]);
				if (inv) ExOp::bytereverse(tmptmps);
				f_out.push_back((unsigned int) tmptmps);
				if ( arraysize  > 1) {
					tmptmps = *(unsigned short*)&(curflaglist[flagindex*12+10]);
					if (inv) ExOp::bytereverse(tmptmps);
					f_out.push_back((unsigned int) tmptmps);
				}
				break;
			case 0x00000102:
				f_out.push_back(*(unsigned char*)&(curflaglist[flagindex*12+8]));
				if ( arraysize > 1) f_out.push_back(*(unsigned char*)&(curflaglist[flagindex*12+9]));
				if ( arraysize  > 2) f_out.push_back(*(unsigned char*)&(curflaglist[flagindex*12+10]));
				if ( arraysize  > 3) f_out.push_back(*(unsigned char*)&(curflaglist[flagindex*12+11]));
				break;
		}
		return(f_out);
	}

template< >
	vector<long int> TiffFile::getValue< vector<long int> >(int flagindex){
        int rout;
		int type = flagType(flagindex);
		vector<long int> f_out;
		unsigned int tmptmp;
		unsigned short tmptmps;
		unsigned int i = *(unsigned int*)&(curflaglist[flagindex*12+8]);
		unsigned int arraysize = (*(unsigned int*)&(curflaglist[flagindex*12+4]));
		if (inv) {ExOp::bytereverse(i);ExOp::bytereverse(arraysize);}
		switch(type){
			case 0x00000004:
				fseek(f, i, SEEK_SET);
				for(i=0;i< arraysize ;i++) {
					if ((rout = fread(&tmptmp,sizeof(unsigned int), 1, f)) != 1) return f_out;
					if (inv) ExOp::bytereverse(tmptmp);
					f_out.push_back((long int)tmptmp);
				}
				break;
			case 0x00000003:
				fseek(f, i, SEEK_SET);
				for(i=0;i< arraysize ;i++) {
					if ((rout = fread(&tmptmps,sizeof(unsigned short), 1, f)) != 1) return f_out;
					if (inv) ExOp::bytereverse(tmptmps);
					f_out.push_back((long int)tmptmps);
				}
				break;
			case 0x00000104:
				f_out.push_back((long int)i);
				break;
			case 0x00000103:
				tmptmps = *(unsigned short*)&(curflaglist[flagindex*12+8]);
				if (inv) ExOp::bytereverse(tmptmps);
				f_out.push_back((long int) tmptmps);
				if ( arraysize  > 1) {
					tmptmps = *(unsigned short*)&(curflaglist[flagindex*12+10]);
					if (inv) ExOp::bytereverse(tmptmps);
					f_out.push_back((long int) tmptmps);
				}
				break;
			case 0x00000102:
				f_out.push_back((long int)*(unsigned char*)&(curflaglist[flagindex*12+8]));
				if ( arraysize > 1) f_out.push_back((long int)*(unsigned char*)&(curflaglist[flagindex*12+9]));
				if ( arraysize  > 2) f_out.push_back((long int)*(unsigned char*)&(curflaglist[flagindex*12+10]));
				if ( arraysize  > 3) f_out.push_back((long int)*(unsigned char*)&(curflaglist[flagindex*12+11]));
				break;
		}
		return(f_out);
	}


	void TiffFile::addopt_Description(Vector< char* > & opt, char const * const descript) const{
		unsigned int j = opt.getSize();
		unsigned int i = strlen(descript);
		opt.push_back(new char[i+5]);
		memcpy(opt[j], "DESC", sizeof(char)*4);
		memcpy(opt[j]+4, descript, sizeof(char)*(i+1));
	}

	void TiffFile::addopt_Xscale(Vector< char* > & opt, double min_x, double max_x) const{
		unsigned int j = opt.getSize();
		opt.push_back(new char[4 + sizeof(double)*2]);
		memcpy(opt[j], "RNGX", sizeof(char)*4);
		memcpy(opt[j]+4, &min_x, sizeof(char)*4);
		memcpy(opt[j]+4+sizeof(double),  &max_x, sizeof(char)*4);
	}
	void TiffFile::addopt_Yscale(Vector< char* > & opt, double min_y, double max_y) const{
		unsigned int j = opt.getSize();
		opt.push_back(new char[4 + sizeof(double)*2]);
		memcpy(opt[j], "RNGY", sizeof(char)*4);
		memcpy(opt[j]+4, &min_y, sizeof(char)*4);
		memcpy(opt[j]+4+sizeof(double),  &max_y, sizeof(char)*4);
	}

    void TiffFile::startFrameWrite(){
        current = new WriteScope();
    }

    void TiffFile::endFrameWrite(){
        delete[](current);

    }
    void TiffFile::savePLYasTIFF(const char * path){
        FILE* f = fopen(path, "r+");
        char buffer[65536];
        unsigned int nbvertex;
        unsigned int nbface;
        // parse header
        while(1 == fscanf(f," %s",buffer)){
            if (strcmp(buffer, "element") == 0){
                fscanf(f," %s",buffer);
                if (strcmp(buffer, "vertex") == 0) fscanf(f," %i", &nbvertex);
                else if (strcmp(buffer, "face") == 0)fscanf(f," %i", &nbface);
            } else if (strcmp(buffer, "end_header") == 0) break;
            }

        Tuple<unsigned int, 2> vertex_coor;
        vertex_coor[0] = 8;
        vertex_coor[1] = nbvertex;
        DataGrid<float, 2> vertexmap(vertex_coor);
        unsigned int i,j;
        float* buf;
        float maxcoor =0.0f;
        for(i=0;i<nbvertex;i++) {
            buf = vertexmap.data + i * 8;
            fscanf(f," %f %f %f %f %f %f %f %f", buf+1, buf+2, buf, buf+4, buf+5, buf+3, buf+6, buf+7);
            for(j=0;j<3;j++) if (maxcoor < fabs(buf[j])) maxcoor = fabs(buf[j]);
        }

        for(i=0;i<nbvertex;i++) {
            buf = vertexmap.data + i * 8;
            for(j=0;j<3;j++) buf[j] /= maxcoor;
            /*
            if (fabs(buf[4]) >= 0.9f){
                buf[6] = (buf[0] + 1.0f) * ((buf[4]  > 0.0f) ? 0.5f : -0.5f);
                buf[7] = (buf[2] + 1.0f) * 0.5f;
            }else if (fabs(buf[5]) >= 0.999f){
                buf[6] = (buf[0] + 1.0f) * ((buf[5]  > 0.0f) ? 0.5f : -0.5f);
                buf[7] = (buf[1] + 1.0f) *0.5f;
            }else if (fabs(buf[3]) >= 0.999f){
                buf[6] = (buf[2] + 1.0f) * ((buf[3]  > 0.0f) ? 0.5f : -0.5f);
                buf[7] = (buf[1] + 1.0f) *0.5f;
            }else if (fabs(buf[4]) < sqrt(0.095f)){
                buf[6] = ((float)rand()) / RAND_MAX;
                buf[7] = (buf[1] + 1.0f) *0.5f;
            }*/
        }
        this->put(vertexmap, (float) -1.0f, (float) 1.0f);

        Vector<unsigned int> triangles;
        Vector<unsigned int> quads;
        for(i=0;i<nbface;i++) {
            fscanf(f," %i", &(vertex_coor[0]));
            switch(vertex_coor[0]){
                case 3: for(vertex_coor[1]=0;vertex_coor[1]<3;vertex_coor[1]++) {fscanf(f," %i", &(vertex_coor[0])); triangles.push_back(vertex_coor[0]); } break;
                case 4: for(vertex_coor[1]=0;vertex_coor[1]<4;vertex_coor[1]++) {fscanf(f," %i", &(vertex_coor[0])); quads.push_back(vertex_coor[0]); } break;
                default: fscanf(f,"%*[^\n]"); // flush line
            }
        }
        DataGrid<unsigned int, 2> index_map;
        if (triangles.getSize() != 0){
            vertex_coor[0] =3;
            vertex_coor[1] =triangles.getSize() / 3;
            index_map.setSizes(vertex_coor);
            memcpy(index_map.data, triangles.darray, triangles.getSize() * sizeof(unsigned int) );
            this->put(index_map, (unsigned int) 0, (unsigned int) nbvertex-1);
        }
        if (quads.getSize() != 0){
            vertex_coor[0]= 4;
            vertex_coor[1]= quads.getSize() / 4;
            index_map.setSizes(vertex_coor);
            memcpy(index_map.data, quads.darray, quads.getSize() * sizeof(unsigned int) );
            this->put(index_map, (unsigned int) 0, (unsigned int) nbvertex-1);
        }
        fclose(f);
        printf("Conversion done!, %i vertice, %i triangles %i quads\n", nbvertex, triangles.getSize()/3, quads.getSize()/4);
    }



	void Bluestein::setSize(unsigned int n_size){
		tsize= n_size;
		pre_pow2 =0;
		if (n_size == 0) {static_warning_handdle << LFH_WARNING_UNEXPECTED_INPUT; return;}
		while((tsize & 1) == 0) {pre_pow2++;tsize >>= 1;}
		if (tsize == 1){ // pure power of 2!
			post_pow2 =0;
			mult = 0;
		}else{
			mult = tsize;
			unsigned int t = (tsize << 1) + 1;
			for(post_pow2=0;t != 0;post_pow2++) t >>=1;


			post_pow2 += pre_pow2;
			mult <<= pre_pow2;
			pre_pow2 =0;

			bluewindow.setSize(1 << post_pow2);

			unsigned int i;

			for(i=0; i < mult;i++) {
				double ang = i*i*M_PI/ mult;
				bluewindow[i] =  mycomplex(cos(ang),sin(ang));
			}
			for(;((i + mult-1)  >> post_pow2) == 0;i++) bluewindow[i] =  ExCo< mycomplex >::mkZero();

			for(;(i >> post_pow2) == 0;i++) {
				double ang = ((1 << post_pow2)-i)*((1 << post_pow2)-i)*M_PI/ mult;
				bluewindow[i] =  mycomplex(cos(ang),sin(ang));
			}
			pow2_FFT_routine(bluewindow.darray, post_pow2);
		}
		tsize= n_size;
	}

	unsigned int Bluestein::getBufferSize()const{
		if (mult == 0) return tsize;
		else return tsize - mult + (1 << post_pow2);
	}

	void Bluestein::show(FILE* f, int level)const{
		fprintf(f, "Bluestein! %i,%i,%i,%i;\n", pre_pow2, mult, post_pow2, getBufferSize());
	}

	FoldedGaussianDistribution::FoldedGaussianDistribution() {}
	FoldedGaussianDistribution::FoldedGaussianDistribution(double _mean, double _std): mean(_mean),i_std(sqrt(0.5f)/_std){}

	void FoldedGaussianDistribution::operator()(double &_out_prob, double &obs) const{

		double tmp = (mean - obs) * i_std;
		double tmp2 = (mean + obs) * i_std;
		_out_prob = i_std * (exp(-tmp * tmp) + exp(-tmp2*tmp2)) * (1.0f / sqrt(M_PI));

	}


	double FoldedGaussianDistribution::LL(const double &obs) const{
		double tmp = (mean - obs) * i_std;
		double tmp2 = (mean + obs) * i_std;
		return log(i_std / sqrt(M_PI)) * log(exp(-tmp * tmp) + exp(-tmp2*tmp2));
		// X = log(COSH X + SINH X)

	}



	void FoldedGaussianDistribution::EMinit(){
	}

	void FoldedGaussianDistribution::EMAlphainit(double){
	}

	void FoldedGaussianDistribution::EMregist( const double &instance, const double prob){
	}

	double FoldedGaussianDistribution::EMfinit(){
		return(0.0f);
	}

	FoldedGaussianDistribution FoldedGaussianDistribution::EMnext_exec(double alpha) const{
		FoldedGaussianDistribution ha;
		return ha;
	}

	void FoldedGaussianDistribution::EMclear(){
	}

	void FoldedGaussianDistribution::show(FILE* o) const{
	}

void VarGammaDistribution::operator()(double &_out_prob, double &obs) const{
    _out_prob =exp( this->LL(obs));
}

double VarGammaDistribution::LL(const double&)const{
    return 0.0;
}


	void GaussianDistributionV::setSizes(unsigned int s){
		if (s == nbsizes) return;
		if (nbsizes != 0){
			delete[](mean);
		}
		mean = new double[s*2]; // mean and var eigevals
		unsigned int coor[2];
		coor[0] =s;
		coor[1] =s;
		eigenvec.setSizes(coor);
		if (stat_Scope != NULL) {
			delete[](stat_Scope);
			stat_Scope = NULL;
		}
		nbsizes = s;
	}


	void GaussianDistributionV::operator()(double & prob, double * & dd) const{

	}

	double GaussianDistributionV::LL( double * const & dd) const{
		return(0.0f);
	}
	void GaussianDistributionV::EMinit(){
		if (stat_Scope == NULL) stat_Scope = new double[2 +  ((nbsizes*(nbsizes+3))/2)];
		for(unsigned int i=0;i<2 + ((nbsizes*(nbsizes+3))/2) ;i++) stat_Scope[i] =0.0f;
	}
	void GaussianDistributionV::EMAlphainit(double){
		if (stat_Scope == NULL) stat_Scope = new double[2 + ((nbsizes*(nbsizes+3))/2)];
		for(unsigned int i=0;i<2 + ((nbsizes*(nbsizes+3))/2) ;i++) stat_Scope[i] =0.0f;
	}
	void GaussianDistributionV::EMregist( double * const  &instance, const double prob){
		unsigned int i,j;

		for(i=0;i<nbsizes;i++) if (!ExOp::isValid(instance[i])) break;
		if (i<nbsizes) return;

		stat_Scope[0] += prob;
		stat_Scope[1] += prob*prob;

		for(i=0;i<nbsizes;i++) stat_Scope[2+i] += instance[i] * prob;
		unsigned int k = nbsizes+2;
		for(i=0;i<nbsizes;i++){
			for(j=i;j<nbsizes;j++,k++) stat_Scope[k] += instance[i] * instance[j] * prob;
		}
	}


	double GaussianDistributionV::EMfinit(){
		unsigned int i;
		unsigned int coor[2];
		unsigned int acoor[2];
	//	for(i=0;i< 2+((nbsizes*(nbsizes+3))/2) ;i++) printf("%f\t", stat_Scope[i]);
		double num = (stat_Scope[1] - stat_Scope[0]*stat_Scope[0]) /2;

		double tmp = 1.0f / stat_Scope[0];
		for(i=0;i<nbsizes;i++) {
			mean[i] = stat_Scope[i+2] * tmp;
		}

	//	DataGrid<double, 2> mat;
		coor[0] = nbsizes;
		coor[1] = nbsizes;
		eigenvec.setSizes(coor);

		i = 2 + nbsizes;
		for(coor[0]=0;coor[0]<nbsizes;coor[0]++){
		//	ExOp::show(coor);
			tmp = stat_Scope[0] * stat_Scope[i] - stat_Scope[2+coor[0]] * stat_Scope[2+coor[0]];
			if (tmp <= 0) tmp = 0.0f;
			coor[1]=coor[0];
			eigenvec(coor) = tmp;
			acoor[1] =coor[0];
			for(coor[1]++,i++;coor[1]<nbsizes;coor[1]++, i++){
		//		ExOp::show(coor);
				tmp = stat_Scope[0] * stat_Scope[i] - stat_Scope[2+coor[1]] * stat_Scope[2+coor[0]];
				eigenvec(coor) = tmp;
				acoor[0] =coor[1];
		//		ExOp::show(acoor);
				eigenvec(acoor) = tmp;
			}
		}
//		(mat / -2.0f * num) is the covariance Matrix!!!

		eigenvec *= (-0.5f /  num);

	//	Vector<double> eig;
	//	eigenvec.makeDiagonalizer_ofinverse(mat, eig);

	//	ExOp::show(mat);
	//	for(i=0;i<nbsizes;i++) {
	//		mean[i+nbsizes] = eig[i] * num;
	//	}

		return (0.0f); // todo
	}

	void GaussianDistributionV::EMclear(){
		if (stat_Scope != NULL) {
			delete[](stat_Scope);
			stat_Scope = NULL;
		}
	}

	void GaussianDistributionV::show(FILE* o) const{
		fprintf(o,"mean: ");
		unsigned int i;
		for(i=0;i<nbsizes;i++) fprintf(o,"%f%c", mean[i], i +1 == nbsizes ? '\n' : '\t');
		fprintf(o,"EigenVar: ");
		for(i=0;i<nbsizes;i++) fprintf(o,"%f%c", mean[i+nbsizes], i +1 == nbsizes ? '\n' : '\t');


		ExOp::show(eigenvec,o,0);

	}


	void GaussianDistributionV::save(FILE* f) const{
		//		fwrite(&normfactor,sizeof(double),1,f);
		//		fwrite(&mean,sizeof(Tuple<double, nbchannels>),1,f);
		//		fwrite(&ihvar,sizeof(Matrix<double, nbchannels, nbchannels>),1,f);
	}

	void GaussianDistributionV::load(FILE* f, unsigned int size){
		//		fread(&normfactor,sizeof(double),1,f);
		//		fread(&mean,sizeof(Tuple<double, nbchannels>),1,f);
		//		fread(&ihvar,sizeof(Matrix<double, nbchannels, nbchannels>),1,f);
	}

	GradientSearchScope::GradientSearchScope(const GradientSearchScope& in) : size(in.size),lastvalue(in.lastvalue), expdiff(in.expdiff), log_alpha(in.log_alpha){
		if (in.lastderiv) {lastderiv = new double[size]; memcpy(lastderiv,in.lastderiv,sizeof(double)*size);}
		else lastderiv = NULL;
	}

	GradientSearchScope& GradientSearchScope::operator=(const GradientSearchScope& in){
		size = in.size;
		if (in.lastderiv) {lastderiv = new double[size]; memcpy(lastderiv,in.lastderiv,sizeof(double)*size);}
		else lastderiv = NULL;
		lastvalue = in.lastvalue;
		expdiff = in.expdiff;
		log_alpha = in.log_alpha;
		return(*this);
	}

	double GradientSearchScope::operator()()const{return(exp(log_alpha));}
	void GradientSearchScope::punish(double log_mag){log_alpha -= log_mag; expdiff *= exp(-log_mag);}

	void GradientSearchScope::init(double initalpha, int in_size){
		size = in_size;
		if (lastderiv != NULL) delete[](lastderiv);
		lastderiv = NULL;
		log_alpha = log(initalpha);
	}

	void GradientSearchScope::registAscent(const double &value, double const * const deriv, double mod_fact){
		int i;
		if (lastderiv == NULL){ // first time, dont update alpha
			lastderiv = new double[size];
			lastvalue = value;
			expdiff  = 0;
			memcpy(lastderiv,deriv,sizeof(double)*size);
			for(i=0;i<size;i++) expdiff  += deriv[i]*deriv[i];
			expdiff = lastvalue + exp(log_alpha) * expdiff;
		}else{
			// update alpha!
			double sum, suma ,sumb, tmp, fact;
			tmp = (value - expdiff); sum = (tmp * tmp);
			tmp = (lastvalue - expdiff); suma = (tmp * tmp); // magnitude of prediction
			tmp = (value - lastvalue); sumb = tmp * tmp;
			fact = exp(- 0.7f * sum * pow(suma*sumb,-0.5f));
			sum =0;
			suma =0;
			sumb =0;
			lastvalue = value;
			expdiff = 0;
			for(i=0;i<size;i++) {tmp = lastderiv[i] - deriv[i]; sum += tmp*tmp; suma += deriv[i]*deriv[i]; sumb += lastderiv[i] * lastderiv[i]; expdiff += deriv[i]*deriv[i];}
			log_alpha += fact + exp(- 0.7f * sum * pow(suma*sumb,-0.5f)) - 1.0f ; // multiply by exp(1) or exp(-1) in extreme cases
			expdiff = lastvalue + exp(log_alpha) *expdiff;
		}
	}
	void GradientSearchScope::registDescent(const double &value, double const * const deriv, double mod_fact){
		int i;
		if (lastderiv == NULL){ // first time, dont update alpha
			lastderiv = new double[size];
			expdiff = 0;
			for(i=0;i<size;i++) expdiff += deriv[i]*deriv[i];
		}else{
			double sump, modif;

			sump =0;
			expdiff = 0;
			for(i=0;i<size;i++) {expdiff += deriv[i] * deriv[i]; sump += deriv[i] * lastderiv[i];}

	//		printf("F[x] = %e; dF[x]/dx = %e; dF[x]/dy = %e; F[x+s]-F[x] = %e; GR = %e; old deriv = %e; new deriv = %e \n", lastvalue, lastderiv[0], lastderiv[1], value - lastvalue, -lexpdiff * exp(log_alpha) * mod_fact, lexpdiff , sump);
	//		printf("curvature! = %e, stepsize %e ,log alpha incr : %e\n", 1.0f / fabs(lexpdiff - sump), exp(log_alpha) * (sqrt(lexpdiff) * mod_fact), 0.5f * log(lexpdiff) + log(mod_fact) - log(fabs(lexpdiff - sump)) );


			//		modif = -0.5f * log( fabs( 1.0f - (sump / expdiff) )); // geometric mean update
					modif = -0.5f * log( 0.1f + fabs( 1.0f - (sump / expdiff) )); // geometric mean update

			if (ExCo<double>::isValid(modif)){
				log_alpha +=  modif;
			}
	//		printf("curvature! = %e   %e\n",log_alpha, modif);

			// F[-S/2] = lastvalue, F[S/2] = value F'[-S] = lexpdiff  F[S] = expdiff

			// F[-S/2] = lastvalue, F[S/2] = value F'[-S] = lexpdiff  F[S] = expdiff

			//  0.5f * (lastvalue + value) = A + C * S^2 // useless!
			//  0.5f * (-lastvalue + value) = B*S + D * S^3
			//  0.5f * (lexpdiff + expdiff) = B * 3D*S^2
			//   (-lexpdiff + expdiff) = C*S //!

			// 1/2C ~=  2.0f / (-lexpdiff + expdiff)




		}
		lastvalue = value;
		memcpy(lastderiv,deriv,sizeof(double)*size);
	}

	double GradientSearchScope::updateAscent(const double &value, double * guess,  double const * const deriv){// returns a safe bound on the parameter changes;
		int i;
		double sums,  sump, modif;
		if (lastderiv == NULL){ // first time, dont update alpha
			lastderiv = new double[size];
			expdiff = 0.0f;
			for(i=0;i<size;i++) {expdiff += deriv[i] * deriv[i];}
		}else if (value <= lastvalue){
               // reject point!, use cubic fit to estimate minimum between previous point an this point F[0] = lastvalue F[S] = value F'[0]

                sump =0;
                for(i=0;i<size;i++) {sump += deriv[i] * lastderiv[i];}
         //       printf("super scope! f0 = %f f1 = %f d0 = %f d1 = %f\n", lastvalue, value, expdiff, sump);
                // distance between points:   alpha *expdiff

                double tmp[3];
                sums = exp(log_alpha);
                tmp[0] = -lastvalue -  sums* expdiff + value; //objective
                tmp[1] = sump - expdiff; //objective // sqrt(expdiff) premulti
          //      printf("%f,%f\n", tmp[0],tmp[1]);

                tmp[2] = -(3.0f * tmp[0] - tmp[1]* sums);
                tmp[0] = 3.0f * (2.0f * tmp[0]- tmp[1]* sums); // sqrt(expdiff) premulti
           //     printf("x2 =%f  x3 = %f in absolute direction\n", tmp[2], tmp[3]);

          //      printf("deriv discriminant, be positive :/ = %f\n", tmp[4]);
                tmp[1] = -tmp[2] / tmp[0];
                if (fabs(tmp[1]) > 1000000000.0f){ // pivot is too large!
                    tmp[0] = sums * expdiff / (tmp[2] * 2.0f);

                }else{
                tmp[0] = tmp[1] + (sqrt(tmp[2] * tmp[2] + tmp[0] * expdiff * sums) / tmp[0]);
           //     printf("altroot: %f\n", tmp[5] - (sqrt(tmp[4]) / tmp[3]));
                }
                printf("in [0,1] range, desired root = %f\n", tmp[0]); // pick same sign as x3

          //      if ((tmp[0] < 0.0f)||(tmp[0]> 1.0f)) printf("NO way!!!!\n");

                //printf("recover guess: ");
                sums *= (1.0f - tmp[0]);
                for(i=0;i<size;i++) guess[i] -= sums * lastderiv[i];
                // alpha = alpha * tmp8 -> logalpha += log(tmp[8])
                log_alpha += log(tmp[0]);
                return(ExCo<double>::mkMaximum());
        }else{
			sums = expdiff;
			sump =0;
			expdiff = 0;
			for(i=0;i<size;i++) {expdiff += deriv[i] * deriv[i]; sump += deriv[i] * lastderiv[i];}
			//		modif = -0.5f * log( fabs( 1.0f - (sump / sums) )); // geometric mean update
			modif = -0.5f * log( 0.1f + fabs( 1.0f - (sump / sums) )); // geometric mean update
			if (ExCo<double>::isValid(modif)){
				log_alpha +=  modif;
			}
		}
		lastvalue = value;
		memcpy(lastderiv,deriv,sizeof(double)*size);
		sums = exp(log_alpha);
		modif = 0.0f;
		for(i=0;i<size;i++) {sump = guess[i] ; guess[i] += sums * deriv[i]; modif += (fabs(guess[i]) > 1.2e-200) ? fabs((guess[i] - sump)/ guess[i]) : 0.0f ;}
		return(modif /size);
	}



	double GradientSearchScope::updateDescent(const double &value, double * guess , double const * const deriv){ // returns a safe bound on the parameter changes;
		int i;
		double sums,  sump, modif;
		if (lastderiv == NULL){ // first time, dont update alpha
			lastderiv = new double[size];
			expdiff = 0.0f;
			for(i=0;i<size;i++) {expdiff += deriv[i] * deriv[i];}
		}else if (value >= lastvalue) {
               // reject point!, use cubic fit to estimate minimum between previous point an this point F[0] = lastvalue F[S] = value F'[0]
                sump =0;
                for(i=0;i<size;i++) {sump += deriv[i] * lastderiv[i];}
         //       printf("super scope! f0 = %f f1 = %f d0 = %f d1 = %f\n", lastvalue, value, expdiff, sump);
                // distance between points:   alpha *expdiff

                double tmp[3];
                sums = exp(log_alpha);
                tmp[0] = lastvalue -  sums* expdiff - value; //objective
                tmp[1] = sump - expdiff; //objective // sqrt(expdiff) premulti

                tmp[2] = -(3.0f * tmp[0] - tmp[1]* sums);
                tmp[0] = 3.0f * (2.0f * tmp[0]- tmp[1]* sums); // sqrt(expdiff) premulti
           //     printf("x2 =%f  x3 = %f in absolute direction\n", tmp[2], tmp[3]);

          //      printf("deriv discriminant, be positive :/ = %f\n", tmp[4]);
                tmp[1] = -tmp[2] / tmp[0];
                if (fabs(tmp[1]) > 1000000000.0f){ // pivot is too large!
                    tmp[0] = sums * expdiff / (tmp[2] * 2.0f);
                }else{
                tmp[0] = tmp[1] + (sqrt(tmp[2] * tmp[2] + tmp[0] * expdiff * sums) / tmp[0]);
           //     printf("altroot: %f\n", tmp[5] - (sqrt(tmp[4]) / tmp[3]));
                }
               // printf("in [0,1] range, desired root = %f\n", tmp[0]); // pick same sign as x3

                if ((tmp[0] < 0.0f)||(tmp[0]> 1.0f)) printf("NO way!!!!\n");

                //printf("recover guess: ");
                sums *= (1.0f - tmp[0]);
                for(i=0;i<size;i++) guess[i] += sums * lastderiv[i];
                // alpha = alpha * tmp8 -> logalpha += log(tmp[8])
                log_alpha += log(tmp[0]);
                return(ExCo<double>::mkMaximum());
        }else{
			sums = expdiff;
			sump =0;
			expdiff = 0;
			for(i=0;i<size;i++) {expdiff += deriv[i] * deriv[i]; sump += deriv[i] * lastderiv[i];}
			//		modif = -0.5f * log( fabs( 1.0f - (sump / sums) )); // geometric mean update
	//		printf("%e expected, %e obtained\n", );
			modif = -0.5f * exp(log_alpha) - ((value - lastvalue) / (sump + sums)); // cubic error
			modif = -0.5f * log( fabs( modif) + fabs( 1.0f - (sump / sums) )); // geometric mean update
			if (ExCo<double>::isValid(modif)){
				log_alpha +=  modif;
			}

		}
		lastvalue = value;
		memcpy(lastderiv,deriv,sizeof(double)*size);
		sums = exp(log_alpha);
		modif = 0.0f;
		for(i=0;i<size;i++) {sump = guess[i] ; guess[i] -= sums * deriv[i]; modif += (fabs(guess[i]) > 1.2e-200) ? fabs((guess[i] - sump)/ guess[i]) : 0.0f ;}
		return(modif /size);
	}
CurvatureSearchScope::CurvatureSearchScope(const CurvatureSearchScope& in) : size(in.size),lastvalue(in.lastvalue), expdiff(in.expdiff){
    if (in.lastderiv) {lastderiv = new double[size*4]; memcpy(lastderiv,in.lastderiv,sizeof(double)*size*4);
        curvature.setSize(size & 1 ? ((size >> 1)+1)*(size) : (size >> 1)*(size+1)); memcpy(lastderiv,in.lastderiv,sizeof(double)*size);
    }
    else lastderiv = NULL;
}
CurvatureSearchScope& CurvatureSearchScope::operator=(const CurvatureSearchScope& in){
    size = in.size;
    if (in.lastderiv) {lastderiv = new double[size*4]; memcpy(lastderiv,in.lastderiv,sizeof(double)*size*4);}
    else lastderiv = NULL;
    lastvalue = in.lastvalue;
    expdiff = in.expdiff;
    alpha = in.alpha;
    shrinkcount = in.shrinkcount;
    curvature = in.curvature;
return(*this);}
void CurvatureSearchScope::init(unsigned int in_size, double initalpha, double _epsilon){
    epsilon = _epsilon;
    size = in_size;
    if (lastderiv != NULL) delete[](lastderiv);
    lastderiv = NULL;
    curvature.setSize(size & 1 ? ((size >> 1)+1)*(size) : (size >> 1)*(size+1));
    unsigned int i,j,k;
    for(i=0,k=0;i<size;i++) {
        curvature[k] = initalpha;
        for(j=i+1,k++;j<size;j++,k++) curvature[k] = 0.0f;
    }
    shrinkcount =0;
}
	void CurvatureSearchScope::init(unsigned int in_size, const Trianglix<double> &dderiv, double _epsilon){
        epsilon = _epsilon;
        size = in_size;
		if (lastderiv != NULL) delete[](lastderiv);
		lastderiv = NULL;
		curvature.setSize(size & 1 ? ((size >> 1)+1)*(size) : (size >> 1)*(size+1));
        unsigned int i,j,k;
        double scale_prior = 0.0;
        uint32_t nbnan =0u;
        Trianglix<double> invdd = dderiv.mkInverse();
        for(i=0,k=0;i<size;i++) {
            curvature[k] = 1.0 / fabs(dderiv[i]);
            if (ExOp::isValid(curvature[k])) scale_prior -= log(fabs(dderiv[i]));
            else nbnan++;
            for(j=i+1,k++;j<size;j++,k++) curvature[k] = 0.0;
        }
        if (nbnan != 0){
            if (nbnan == size){
                printf("Provided an double deriv vector with zeros/nans only...\n");
                exit(1);
            }
            scale_prior = exp(scale_prior / (size - nbnan));


       //    printf("Fixed %i ddervi with a prior... %e\n", nbnan,  scale_prior);
            for(i=0,k=0;i<size;i++) {
                if (!ExOp::isValid(curvature[k])) curvature[k] = scale_prior;
                k += size - i;
            }
        }
        shrinkcount =0;
	}
void CurvatureSearchScope::init(unsigned int in_size, const double* dderiv, double _epsilon){
    epsilon = _epsilon;
    size = in_size;
    if (lastderiv != NULL) delete[](lastderiv);
    lastderiv = NULL;
    curvature.setSize(size & 1 ? ((size >> 1)+1)*(size) : (size >> 1)*(size+1));
    unsigned int i,j,k;
    double scale_prior = 0.0;
    uint32_t nbnan =0u;
    for(i=0,k=0;i<size;i++) {
        curvature[k] = 1.0f / fabs(dderiv[i]); // TMPTMP! checking that it is stable
        if (ExOp::isValid(curvature[k])) scale_prior -= log(fabs(dderiv[i]));
        else nbnan++;
        for(j=i+1,k++;j<size;j++,k++) curvature[k] = 0.0;
    }
    if (nbnan != 0){
        if (nbnan == size){
            printf("Provided an double deriv vector with zeros/nans only...\n");
            exit(1);
        }
        scale_prior = exp(scale_prior / (size - nbnan));
        for(i=0,k=0;i<size;i++) {
            if (!ExOp::isValid(curvature[k])) curvature[k] = scale_prior;
            k += size - i;
        }
    }
    shrinkcount =0;
}
	void CurvatureSearchScope::initMK2(unsigned int in_size, const double* dderiv){ // need to handdle dderiv >= 0...
        size = in_size;

		if (lastderiv != NULL) delete[](lastderiv);
		lastderiv = NULL;

		curv.setSize(size);
        unsigned int i,j,k;
        double scale_prior = 0.0;
        uint32_t nbnan =0;
        for(i=0,k=0;i<size;i++,k++) {
            for(j=0;j<i;j++) curv.data[k++] = 0.0;
            curv.data[k] = 1.0 / fabs(dderiv[i]);
            if (ExOp::isValid(curv.data[k])) scale_prior -= log(fabs(dderiv[i]));
            else nbnan++;
        }
        if (nbnan != 0){
            if (nbnan == size){
                printf("Provided an double deriv vector with zeros/nans only...\n");
                exit(1);
            }
            scale_prior = exp(scale_prior / (size - nbnan));


       //    printf("Fixed %i ddervi with a prior... %e\n", nbnan,  scale_prior);
            for(i=0,k=0;i<size;i++) {
                if (!ExOp::isValid(curv.data[k])) curv.data[k] = scale_prior;
                k += i+2;
            }
        }
        shrinkcount =0;
	}

	void CurvatureSearchScope::updateRandom(double * guess){
        // make a random move, & throw out derivative data..
        uint32_t i;
        if (lastderiv == NULL) myexit("dont use randow twice in a row...");
        double difmag =0.0f;
        for(i=0;i<size;i++) difmag += log(fabs(guess[i] - (lastderiv[size+i])));
        difmag =exp(difmag /size);

        //double step;
        for(i=0;i<size;i++){
            guess[i] += difmag * Pvalue_to_stdnorm( (0.5 + rand()) / RAND_MAX );
        }
        delete[](lastderiv); lastderiv = NULL;
	}

void CurvatureSearchScope::updateGotNAN(double * guess){
    if (lastderiv == NULL) {
        fprintf(stderr, "initial guess has function value NaN, provide a better guess...\n");
        exit(1);
    }else{
        memcpy(guess,lastderiv+size,sizeof(double)*size);
        uint32_t i,j,k;
        for(i=0,k=0;i<size;i++) {
            curvature[k] *= 0.0625;
            guess[i] += curvature[k] * lastderiv[i];//printf("%e\t",curvature[k] );
            for(j=i+1,k++;j<size;j++,k++){
                curvature[k] *= 0.0625;
                guess[j] += curvature[k] * lastderiv[i];//printf("%e\t",curvature[k] );
                guess[i] += curvature[k] * lastderiv[j];
            }
        }
    }
}

/** \brief UpdateAscent: updates
 *
 * \param value
 * \param guess
 * \param deriv
 * \return true if convergence reached
 *
 */
bool CurvatureSearchScope::updateAscent(const double &value, double * guess,  double const * const deriv){// returns a safe bound on the parameter changes;
    unsigned int i,j,k;
    double sums,  sump;
    if (lastderiv == NULL) {
        if (!ExOp::isValid(value)) {fprintf(stderr, "initial guess has function value NaN, provide a better guess...\n"); return true;}
        lastderiv = new double[size*4];
        lastderiv[size * 3] = ExCo<double>::mkMaximum();
    }else if (!ExOp::isValid(value)){ // not-a-number eh?
        /*sump = 0.0f;
        for(i=0;i<size;i++) {
            sump += lastderiv[i+ size *2] * lastderiv[i];
        }
        sump = -0.75 / sump; // factor

        for(i=0,k=0;i<size;i++) {
            curvature[k] += sump * lastderiv[i+ size *2] *  lastderiv[i+ size *2];// printf("%f\t", curvature[k]);
            for(j=i+1,k++;j<size;j++,k++){
                curvature[k] += sump * lastderiv[i+ size *2] * lastderiv[j+ size *2]; //printf("%f\t", curvature[k]);
            }
        }*/
       // memcpy(lastderiv + size * 3 ,guess ,sizeof(double)*size);
        //for(i=0;i<size;i++){printf("%e\t", lastderiv[i + size*2]);} printf(" before\n");


        memset(lastderiv + size * 2 ,'\0', sizeof(double)*size);
        sump = 0.0f; lastderiv[size * 3] =0.0;
        for(i=0,k=0;i<size;i++) {
            curvature[k] *= 0.0625f;
            lastderiv[i+ size *2] += curvature[k] * lastderiv[i];//printf("%e\t",curvature[k] );
            for(j=i+1,k++;j<size;j++,k++){
                curvature[k] *= 0.0625f;
                lastderiv[j+ size *2] += curvature[k] * lastderiv[i];
                lastderiv[i+ size *2] += curvature[k] * lastderiv[j];//printf("%e\t",curvature[k] );
            }
            sump += lastderiv[i+ size *2] * lastderiv[i+ size *2];
            guess[i] = lastderiv[i+ size *2] + lastderiv[i + size];
            if (guess[i] != lastderiv[i + size]){
                sums = fabs(lastderiv[i+ size *2] = guess[i] - lastderiv[i + size]) / (fabs(guess[i]) + fabs(lastderiv[i + size]));
                if (lastderiv[size * 3] < sums) lastderiv[size * 3] = sums;
            }else lastderiv[i+ size *2] = 0;
        }
        /*
        for(i=0;i<size;i++){printf("%e\t", lastderiv[size+i]);}
        printf(", norm=%e (safe?) \n", sump);

        for(i=0;i<size;i++){printf("%e\t", guess[i] - lastderiv[size+i]);}
        printf(", norm=%e (fc %e)\n", sump, sump / sums);
        for(i=0;i<size;i++){printf("%e\t", (guess[i] - lastderiv[i + size]) / lastderiv[i + size *3]);}
        */
        //for(i=0;i<size;i++){printf("%e\t", lastderiv[i + size*2]);} printf(" after\n");
        //printf("%e\t%e\tOut of bounds Step %e\n", lastvalue, lastderiv[size * 3], sqrt(sump));
        return false;
    }else{
        // check Wolfe Conditions! so that the hessian estimate remains positive definite

        // value - lastvalue = 0.0001(guess - oldguess) \cdot oldderiv
        // (guess - oldguess) \cdot deriv = 0.9 (guess - oldguess) \cdot oldderiv
        sums =0; sump =0;
        /*for(i =0;i<size;i++){
            printf("%e vs %e\n", guess[i] -  lastderiv[i+ size], lastderiv[i+ size *2]);
        }*/

        // is a maximization; hence, flip signs
        for(i=0;i<size;i++) {sums -= lastderiv[i] * lastderiv[i+ size*2]; // rounding errors would be accounted here!
                                 sump -= deriv[i] * lastderiv[i+ size*2];} //* (guess[i] -  lastderiv[i+ size]
        // these are defined for *minimization* so flip sign!

        lastderiv[size * 3] = fabs( (lastvalue == 0.0) ? value : lastvalue);
        if (lastderiv[size * 3] != 0.0) lastderiv[size * 3] = (fabs((value - lastvalue) / lastderiv[size * 3]) <= pow(0.5, 53)) ? 0.0 : value - lastvalue;

        //printf("wolfe 1: %e < %e  %c\n",  -lastderiv[size * 3], -0.0001f * sums, (lastderiv[size * 3] >= 0.0001f * sums) ? 'Y' : 'N' );
        //printf("wolfe 2: %e < %e  %c\n", -sump, -0.9f * sums , (sump >= 0.9f * sums) ? 'Y' : 'N' );

        if (lastderiv[size * 3] < 0.0001f * sums){  // step too long
            // scaling inverse hessian by 0.25
            // sums is the denominator!

            if (sump < 0.9f * sums) { // step is too short AND too long... so this is it.
                if (value < lastvalue) {
                    //printf("Confused Bad Exit!\n");
                    lastvalue = value;
                    memcpy(guess,lastderiv + size ,sizeof(double)*size);
                }//else printf("Confused Good Exit!\n");
                return true; // not time expand, rounding errors makes derivative crazy I say
            }

            /*memset(lastderiv + size*3 ,'\0', sizeof(double)*size);
            sump=0;
            for(i=0,k=0;i<size;i++) {
                lastderiv[i + size*3]+= curvature[k] * lastderiv[i];//printf("%e\t",curvature[k] );
                for(j=i+1,k++;j<size;j++,k++){
                    lastderiv[j + size*3] += curvature[k] * lastderiv[i];
                    lastderiv[i + size*3] += curvature[k] * lastderiv[j];//printf("%e\t",curvature[k] );
                }
                sump += lastderiv[i + size*3] * lastderiv[i];
            }
            printf("%e vs %e\n", sump, sums);*/
            sums = 0.75 / sums; // the sign of sums is flipped...


            for(i=0,k=0;i<size;i++) {
                curvature[k++] += sums * lastderiv[i+ size *2] * lastderiv[i+ size *2];
                for(j=i+1;j<size;j++){
                    curvature[k++] += sums * lastderiv[i+ size *2] * lastderiv[j+ size *2];
                }
            }
            //for(i=0;i<size;i++){printf("%e\t", lastderiv[i + size*2]);} printf("\n");
            lastderiv[size * 3] =0.0;
            //memcpy(lastderiv + size * 3 ,lastderiv + size * 2 ,sizeof(double)*size);
            memset(lastderiv + size * 2 ,'\0', sizeof(double)*size);

            sump =0.0;
            for(i=0,k=0;i<size;i++) {
                lastderiv[i+ size *2] += curvature[k] * lastderiv[i];//printf("%e\t",curvature[k] );
                for(j=i+1,k++;j<size;j++,k++){
                    lastderiv[j+ size *2] += curvature[k] * lastderiv[i];
                    lastderiv[i+ size *2] += curvature[k] * lastderiv[j];//printf("%e\t",curvature[k] );
                }
                sump += lastderiv[i+ size *2] * lastderiv[i+ size *2];
                guess[i] = lastderiv[i+ size] + lastderiv[i+ size *2];
                if (guess[i] != lastderiv[i + size]){
                    sums = fabs(lastderiv[i+ size *2] = guess[i] - lastderiv[i + size]) / (fabs(guess[i]) + fabs(lastderiv[i + size]));
                    if (lastderiv[size * 3] < sums) lastderiv[size * 3] = sums;
                }else lastderiv[i+ size *2] = 0;
            }
            //for(i=0;i<size;i++){printf("%e\t", lastderiv[i + size*2]);} printf("\n");
            //for(i=0;i<size;i++){printf("%e\t", lastderiv[i + size*2] / lastderiv[i + size*3]);} printf(" fc\n");
            //printf("%e\t%e\tRetract Step %e\n", lastvalue, lastderiv[size * 3], sqrt(sump));
            return false;
        }else if (sump < 0.9f * sums){ // step too short, accept change but not the normal hessian update
            // scaling inverse hessian by 4
            // sums is the denominator!
            //printf("Expand! value change %e\n", value - lastvalue);
            //printf("Expand! %e < 0.9 * %e\n", sump, sums);
            //if (value < lastvalue) printf("Going bad! Cheating! Expand\n");
            if (lastderiv[size * 3] <= pow(0.5, 45)) {
                //printf("Crazy Exit!\n");
                return true;
            } // not time expand, rounding errors makes derivative crazy I say

            sump = 0.0;
            memset(lastderiv + size ,'\0', sizeof(double)*size);
            for(i=0,k=0;i<size;i++) {
                lastderiv[i + size]+= curvature[k] * deriv[i];//printf("%e\t",curvature[k] );
                for(j=i+1,k++;j<size;j++,k++){
                    lastderiv[j + size] += curvature[k] * deriv[i];
                    lastderiv[i + size] += curvature[k] * deriv[j];//printf("%e\t",curvature[k] );
                }
                sump += lastderiv[i + size] * deriv[i];
            }
            /*for(i=0;i<size;i++){printf("%e\t", lastderiv[i + size]);}
            printf("step into deriv %e, fact %e\n", sump, 1.0 / sump);*/

            sump = 3.0 / sump; // factor
            // multiply hessian: H <- (I - a sd)H(I + a ds)
            for(i=0,k=0;i<size;i++) {
                curvature[k++] += sump * lastderiv[i + size] * lastderiv[i + size];// printf("%f\t", curvature[k]);
                for(j=i+1;j<size;j++){
                    curvature[k++] += sump * lastderiv[i + size] *  lastderiv[j + size]; //printf("%f\t", curvature[k]);
                }
            }

            // recompute with new hessian

            memset(lastderiv + size * 2 ,'\0', sizeof(double)*size);
            lastderiv[size * 3] = 0.0;
            sump =0.0;
            for(i=0,k=0;i<size;i++) {
                lastderiv[i+ size *2] += curvature[k] * deriv[i];//printf("%e\t",curvature[k] );
                for(j=i+1,k++;j<size;j++,k++){
                    lastderiv[j+ size *2] += curvature[k] * deriv[i];
                    lastderiv[i+ size *2] += curvature[k] * deriv[j];//printf("%e\t",curvature[k] );
                }
                sump += lastderiv[i+ size *2] * lastderiv[i+ size *2];
                lastderiv[i + size] = guess[i];
                guess[i] += lastderiv[i+ size *2];
                if (guess[i] != lastderiv[i + size]){
                    sums = fabs(lastderiv[i+ size *2] = guess[i] - lastderiv[i + size]) / (fabs(guess[i]) + fabs(lastderiv[i + size]));
                    if (lastderiv[size * 3] < sums) lastderiv[size * 3] = sums;
                }else lastderiv[i+ size *2] = 0;
            }
            //printf("%e\t%e\tExpand  Step %e \n", value, lastderiv[size * 3], sqrt(sump));
            lastvalue = value; memcpy(lastderiv,deriv,sizeof(double)*size);
            return false;
        }else {
            // if condition1 fails, the step is too long: it is impossible to update the hessian
            // if condition2 fails, the step is too short: it is impossible to update the hessian

     //   printf("AlRight!\n");
        //if (value < lastvalue) printf("Going bad! Cheating normal!\n");
        for(i=0;i<size;i++) {
            // lastderiv[i+2*size] = -lastderiv[i+2*size];
            lastderiv[i+3*size] = deriv[i] - lastderiv[i];
        }

        i=0;k=0;

        lastderiv[i+size] = curvature[k] * lastderiv[i+3*size];
        for(j=i+1,k++;j<size;j++,k++){
            lastderiv[j+size] = curvature[k] * lastderiv[i+3*size];
            lastderiv[i+size] += curvature[k] * lastderiv[j+3*size];
        }
        for(i++;i<size;i++) {
            lastderiv[i+size] += curvature[k] * lastderiv[i+3*size];
            for(j=i+1,k++;j<size;j++,k++){
                lastderiv[j+size] += curvature[k] * lastderiv[i+3*size];
                lastderiv[i+size] += curvature[k] * lastderiv[j+3*size];
            }
        }

        sums -= sump;
        sump = -sums;
        for(i=0;i<size;i++) sump += lastderiv[i+size] * lastderiv[i+3*size];
        sums = -1.0f / sums;
        sump *= sums;
   //     printf("K2 = %f\tK1 = %f\n", sump,sums);
        for(i=0,k=0;i<size;i++) {
            curvature[k] += sums * (lastderiv[i+size*2] * (sump * lastderiv[i+size*2] + 2.0f * lastderiv[i+size]));
            for(j=i+1,k++;j<size;j++,k++){
                curvature[k] += sums *(sump * lastderiv[i+size*2]* lastderiv[j+size*2] + (lastderiv[i+size*2]* lastderiv[j+size] + lastderiv[i+size]* lastderiv[j+size*2]));
            }
        }
        lastderiv[size * 3] = 0;

        }
    }
    // has updated curvature by now
    lastvalue = value; memcpy(lastderiv,deriv,sizeof(double)*size);
    memset(lastderiv + size * 2 ,'\0', sizeof(double)*size);
    sump=0.0f;
    for(i=0,k=0;i<size;i++) {
        lastderiv[i+ size *2] += curvature[k] * deriv[i];
        for(j=i+1,k++;j<size;j++,k++){
            lastderiv[j+ size *2] += curvature[k] * deriv[i];
            lastderiv[i+ size *2] += curvature[k] * deriv[j];
        }
        sump += lastderiv[i+ size *2] * lastderiv[i+ size *2];
        lastderiv[i+ size ] = guess[i];
        guess[i] += lastderiv[i+ size *2];
        if (guess[i] != lastderiv[i + size]){
            sums = fabs(lastderiv[i+ size *2] = guess[i] - lastderiv[i + size]) / (fabs(guess[i]) + fabs(lastderiv[i + size]));
            if (lastderiv[size * 3] < sums) lastderiv[size * 3] = sums;
        }else lastderiv[i+ size *2] =0;
    }
    //printf("%e\t%e\tNormal Step %e\n", value, lastderiv[size * 3], sqrt(sump));
return(lastderiv[size * 3] <= pow(0.5, 53));}
double CurvatureSearchScope::getRelErrorBound() const{
    if (lastderiv != NULL) return lastderiv[size * 3];
    else return ExCo<double>::mkMaximum();
}
double CurvatureSearchScope::updateAscentVerbose(const double &value, double * guess,  double const * const deriv, bool debug){// returns a safe bound on the parameter changes;
    unsigned int i,j,k;
    double sums,  sump, modif;
    double chk[4];
    if (lastderiv == NULL) {
        if (!ExOp::isValid(value)) {fprintf(stderr, "initial guess has function value NaN, provide a better guess...\n"); exit(1);}
        lastderiv = new double[size*4];
        modif = ExCo<double>::mkMaximum();
    }else if (!ExOp::isValid(value)){ // not-a-number eh?
        /*sump = 0.0f;
        sums = 0.0f;
        chk[0] = 0.0f;
        for(i=0;i<size;i++) {
            sump += lastderiv[i+ size *2] * lastderiv[i];
            sums += lastderiv[i+ size *2] * lastderiv[i+ size *2];
            chk[0] += lastderiv[i] * lastderiv[i];
        }
        sump = -0.75 / sump; // factor
        //printf("deriv: "); for(i=0;i<size;i++){printf("%e\t", lastderiv[i]);} printf(", norm=%e and factor\n", lastderiv[i], sump);
        // multiply hessian: H <- (I - a sd)H(I + a ds)


        //for(i=0;i<size;i++){printf("%e\t", guess[i]);} printf(", norm=%e\n", sums);
        // multiply hessian: H <- (I - a sd)H(I + a ds)

        for(i=0,k=0;i<size;i++) {
            curvature[k] += sump * lastderiv[i+ size *2] *  lastderiv[i+ size *2];// printf("%f\t", curvature[k]);
            for(j=i+1,k++;j<size;j++,k++){
                curvature[k] += sump * lastderiv[i+ size *2] * lastderiv[j+ size *2]; //printf("%f\t", curvature[k]);
            }
        }*/
       // memcpy(lastderiv + size * 3 ,guess ,sizeof(double)*size);

        // multiply hessian by 0.0625
        memset(lastderiv + size * 2 ,'\0', sizeof(double)*size);
        sump = 0.0f;
        modif = 0.0f;
        for(i=0,k=0;i<size;i++) {
            curvature[k] *= 0.0625f;
            lastderiv[i+ size *2] += curvature[k] * lastderiv[i];//printf("%e\t",curvature[k] );
            for(j=i+1,k++;j<size;j++,k++){
                curvature[k] *= 0.0625f;
                lastderiv[j+ size *2] += curvature[k] * lastderiv[i];
                lastderiv[i+ size *2] += curvature[k] * lastderiv[j];//printf("%e\t",curvature[k] );
            }
            guess[i] = lastderiv[i+ size *2] + lastderiv[i + size];
            sump += lastderiv[i+ size *2] * lastderiv[i+ size *2];
            if (guess[i] != lastderiv[i + size]){
                sums = fabs(guess[i] - lastderiv[i + size]) / (fabs(guess[i]) + fabs(lastderiv[i + size]));
                if (modif < sums) modif = sums;
            }
        }
        /*
        for(i=0;i<size;i++){printf("%e\t", lastderiv[size+i]);}
        printf(", norm=%e (safe?) \n", sump);

        for(i=0;i<size;i++){printf("%e\t", guess[i] - lastderiv[size+i]);}
        printf(", norm=%e (fc %e)\n", sump, sump / sums);
        for(i=0;i<size;i++){printf("%e\t", (guess[i] - lastderiv[i + size]) / lastderiv[i + size *3]);}
        */
        if (debug) printf("value: %e,\tR. step bound %e\tOut of bounds Step\n", lastvalue, modif);
        return ExCo<double>::mkMaximum();
    }else{
        // check Wolfe Conditions! so that the hessian estimate remains positive definite

        // value - lastvalue = 0.0001(guess - oldguess) \cdot oldderiv
        // (guess - oldguess) \cdot deriv = 0.9 (guess - oldguess) \cdot oldderiv
        sums =0; sump =0;
        /*for(i =0;i<size;i++){
            printf("%e vs %e\n", guess[i] -  lastderiv[i+ size], lastderiv[i+ size *2]);
        }*/

        for(i=0;i<size;i++) {sums -= lastderiv[i] * (guess[i] -  lastderiv[i+ size]); // rounding errors would be accounted here!
                                 sump -= deriv[i] * (guess[i] -  lastderiv[i+ size]);} //* (guess[i] -  lastderiv[i+ size]
        // these are defined for *minimization* so flip sign!

        modif = fabs( (value == 0.0) ? lastvalue : value);
        if (modif != 0.0) modif = (((value - lastvalue) / modif) <= pow(0.5, 53)) ? 0.0 : value - lastvalue;

        printf("wolfe 1: %e < %e  %c\n",  -modif, -0.0001f * sums, (modif >= 0.0001f * sums) ? 'Y' : 'N' );
        printf("wolfe 2: %e < %e  %c\n", -sump, -0.9f * sums , (sump >= 0.9f * sums) ? 'Y' : 'N' );

        //

        if (modif  < 0.0001f * sums){  // step too long
            // scaling inverse hessian by 0.25
            // sums is the denominator!

            if (sump < 0.9f * sums){
                printf("Wolfe Confusion!\n");
            }

            sums = 0.75 / fabs(sums);


            for(i=0,k=0;i<size;i++) {
                curvature[k] += sums * lastderiv[i+ size *2] * lastderiv[i+ size *2];
                for(j=i+1,k++;j<size;j++,k++){
                    curvature[k] += sums * lastderiv[i+ size *2] * lastderiv[j+ size *2];
                }
            }
            memset(lastderiv + size * 2 ,'\0', sizeof(double)*size);
            modif =0.0;
            sump=0;
            for(i=0,k=0;i<size;i++) {
                lastderiv[i+ size *2] += curvature[k] * lastderiv[i];//printf("%e\t",curvature[k] );
                for(j=i+1,k++;j<size;j++,k++){
                    lastderiv[j+ size *2] += curvature[k] * lastderiv[i];
                    lastderiv[i+ size *2] += curvature[k] * lastderiv[j];//printf("%e\t",curvature[k] );
                }
                guess[i] += lastderiv[i+ size] + lastderiv[i+ size *2];
                sump += lastderiv[i+ size *2] * lastderiv[i+ size *2];
                if (guess[i] != lastderiv[i + size]){
                    sums = fabs(guess[i] - lastderiv[i + size]) / (fabs(guess[i]) + fabs(lastderiv[i + size]));
                    if (modif < sums) modif = sums;
                }
            }

            if (debug) printf("value: %e,\tr. step bound %e\tContract Step\n", lastvalue, modif);
            return ExCo<double>::mkMaximum();
        }else if (sump < 0.9f * sums){ // step too short, accept change but not the normal hessian update
            // scaling inverse hessian by 4
            // sums is the denominator!
            //printf("Expand! value change %e\n", value - lastvalue);
            //printf("Expand! %e < 0.9 * %e\n", sump, sums);
            shrinkcount =0;
            /*chk[0] = 0.0f;
            chk[1] = 0.0f;
            chk[3] = 0.0f;
            for(i=0;i<size;i++) {
                 chk[0] += (lastderiv[i+size] - guess[i])* (lastderiv[i+size] -guess[i]);
                 chk[1] += lastderiv[i] * lastderiv[i];
                 chk[3] += deriv[i] * deriv[i];
            }
            printf("last step norm is %e, deriv mags %e and %e \n", chk[0],chk[1],chk[3]);*/

            chk[1] =0.0;
             sums = 0.0;
             memset(lastderiv + size *3,'\0', sizeof(double)* size);
            for(i=0,k=0;i<size;i++) {
                lastderiv[i + size *3] += curvature[k] * lastderiv[i];
                for(j=i+1,k++;j<size;j++,k++){
                    lastderiv[i + size *3] += curvature[k] * lastderiv[j];
                    lastderiv[j + size *3] += curvature[k] * lastderiv[i];
                }
                sums += lastderiv[i + size *3] * lastderiv[i + size *3];
                // printf("%i: %e and %e\n", i, guess[i] - lastderiv[i + size] , lastderiv[i + size *3]);
            }


            sump = 0.0;
            memset(lastderiv + size ,'\0', sizeof(double)*size);
            for(i=0,k=0;i<size;i++) {
                lastderiv[i + size]+= curvature[k] * deriv[i];//printf("%e\t",curvature[k] );
                for(j=i+1,k++;j<size;j++,k++){
                    lastderiv[j + size] += curvature[k] * deriv[i];
                    lastderiv[i + size] += curvature[k] * deriv[j];//printf("%e\t",curvature[k] );
                }
                chk[1] += lastderiv[i + size] * lastderiv[i + size];
                sump += lastderiv[i + size] * deriv[i];
            }
            /*for(i=0;i<size;i++){printf("%e\t", lastderiv[i + size]);}
            printf(" norm = \n",chk[1]);
            printf("Constant hessian norm %e  (%e fold)  was %e before\n", chk[1], chk[1]/ chk[0], sums);
            printf("step into deriv %e, fact %e\n", sump, 1.0 / sump);*/

            sump = 3.0 / sump; // factor
            // multiply hessian: H <- (I - a sd)H(I + a ds)
            for(i=0,k=0;i<size;i++) {
                curvature[k] += sump * lastderiv[i + size] * lastderiv[i + size];// printf("%f\t", curvature[k]);
                for(j=i+1,k++;j<size;j++,k++){
                    curvature[k] += sump * lastderiv[i + size] *  lastderiv[j + size]; //printf("%f\t", curvature[k]);
                }
            }

            // recompute with new hessian
            chk[3] = 0.0;
            memcpy(lastderiv + size * 3 ,lastderiv + size ,sizeof(double)*size);
            memset(lastderiv + size * 2 ,'\0', sizeof(double)*size);
            for(i=0,k=0;i<size;i++) {
                lastderiv[i+ size *2] += curvature[k] * deriv[i];//printf("%e\t",curvature[k] );
                for(j=i+1,k++;j<size;j++,k++){
                    lastderiv[j+ size *2] += curvature[k] * deriv[i];
                    lastderiv[i+ size *2] += curvature[k] * deriv[j];//printf("%e\t",curvature[k] );
                }
                guess[i] = lastderiv[i+ size] + lastderiv[i+ size *2];
                chk[3] += lastderiv[i+ size *2] * lastderiv[i+ size *2];
                if (guess[i] != lastderiv[i + size]){
                    sums = fabs(guess[i] - lastderiv[i + size]) / (fabs(guess[i]) + fabs(lastderiv[i + size]));
                    if (modif < sums) modif = sums;
                }
            }
            if (debug) printf("value: %e,\tr. step bound %e\tExpend Step\n", value, modif);
            lastvalue = value; memcpy(lastderiv,deriv,sizeof(double)*size);
            return ExCo<double>::mkMaximum();
        }else { shrinkcount =0;
            // if condition1 fails, the step is too long: it is impossible to update the hessian
            // if condition2 fails, the step is too short: it is impossible to update the hessian

     //   printf("AlRight!\n");

        for(i=0;i<size;i++) {
            lastderiv[i+size] -= guess[i];
            lastderiv[i+3*size] = deriv[i] - lastderiv[i];
        }

        i=0;k=0;

        lastderiv[i+size*2] = curvature[k] * lastderiv[i+3*size];

        for(j=i+1,k++;j<size;j++,k++){
            lastderiv[j+size*2] = curvature[k] * lastderiv[i+3*size];
            lastderiv[i+size*2] += curvature[k] * lastderiv[j+3*size];
   //         printf("%f\t", curvature[k]);
        }
        for(i++;i<size;i++) {
            lastderiv[i+size*2] += curvature[k] * lastderiv[i+3*size];
            for(j=i+1,k++;j<size;j++,k++){
                lastderiv[j+size*2] += curvature[k] * lastderiv[i+3*size];
                lastderiv[i+size*2] += curvature[k] * lastderiv[j+3*size];
            }
        }

        sums -= sump;
        sump = -sums;
        for(i=0;i<size;i++) sump += lastderiv[i+size*2] * lastderiv[i+3*size];
        sums = -1.0f / sums;
        sump *= sums;
   //     printf("K2 = %f\tK1 = %f\n", sump,sums);
        for(i=0,k=0;i<size;i++) {
            curvature[k] += sums * (lastderiv[i+size] * (sump * lastderiv[i+size] - 2.0f * lastderiv[i+size*2]));
            for(j=i+1,k++;j<size;j++,k++){
                curvature[k] += sums *(sump * lastderiv[i+size]* lastderiv[j+size] - (lastderiv[i+size]* lastderiv[j+size*2] + lastderiv[i+size*2]* lastderiv[j+size]));
            }
        }
        modif = 0;

        }
    }
    // has updated curvature by now
    memcpy(lastderiv + size,guess,sizeof(double)*size);
    lastvalue = value; memcpy(lastderiv,deriv,sizeof(double)*size);
    memset(lastderiv + size * 2 ,'\0', sizeof(double)*size);
    for(i=0,k=0;i<size;i++) {
        lastderiv[i+ size *2] += curvature[k] * deriv[i];
        for(j=i+1,k++;j<size;j++,k++){
            lastderiv[j+ size *2] += curvature[k] * deriv[i];
            lastderiv[i+ size *2] += curvature[k] * deriv[j];
        }
        guess[i] = lastderiv[i+ size ] + lastderiv[i+ size *2];
        if (guess[i] != lastderiv[i + size]){
            sums = fabs(guess[i] - lastderiv[i + size]) / (fabs(guess[i]) + fabs(lastderiv[i + size]));
            if (modif < sums) modif = sums;
        }

    }
    if (debug) printf("value: %e,\tR. step bound %e\n", value, modif);
return(modif);}
// use Broyden–Fletcher–Goldfarb–Shanno algorithm
/*
double CurvatureSearchScope::updateAscentMK2(const double &value, double * guess,  double const * const deriv,  double const * const dbl_deriv){// returns a safe bound on the parameter changes;
    unsigned int i,j,k;
    double sums,  sump, modif;
    if (lastderiv == NULL) {
        if (!ExOp::isValid(value)) {fprintf(stderr, "initial guess has function value NaN, provide a better guess...\n"); exit(1);}
        lastderiv = new double[size*4];
        alpha = 0;
    }else if (!ExOp::isValid(value)){
        memcpy(guess,lastderiv+size,sizeof(double)*size);
        modif = 0.0f;
        alpha += 2;
        for(i=0,k=0;i<size;i++) {
            curvature[k] *= 0.0625f;
            guess[i] += curvature[k] * lastderiv[i];//printf("%e\t",curvature[k] );
            for(j=i+1,k++;j<size;j++,k++){
                curvature[k] *= 0.0625f;
                guess[j] += curvature[k] * lastderiv[i];//printf("%e\t",curvature[k] );
                guess[i] += curvature[k] * lastderiv[j];

            }
            modif += (fabs(guess[i]) > 1.2e-100) ? fabs((guess[i] - lastderiv[i+size])/ guess[i]) : 0.0f ;
        }
        return modif;
    }else{

        // check Wolfe Conditions!

        // value - lastvalue = 0.0001(guess - oldguess) \cdot oldderiv
            // (guess - oldguess) \cdot deriv = 0.9 (guess - oldguess) \cdot oldderiv

         //   printf("wolfe 1: %e <= %e  %c\n",  lastvalue - value, 0.0001f * sums, (lastvalue-value <= 0.0001f * sums) ? 'Y' : 'N' );
       //     printf("wolfe 2: %e >= %e  %c\n", sump, 0.9f * sums , (sump>= 0.9f * sums) ? 'Y' : 'N' );

        sums =0; sump =0;
        for(i=0;i<size;i++) {sums -= lastderiv[i] * (guess[i] - lastderiv[i+size]);
                                 sump -= deriv[i] * (guess[i] - lastderiv[i+size]);}
        if (value - lastvalue  < 0.0001f * sums){  // step too long
            // scaling inverse hessian by 0.25
            // sums is the denominator!

        //	printf("Shrink!\n");
            shrinkcount++;


            alpha++;

            sums = -0.75f / fabs(sums);
            sump=0;



            for(i=0;i<size;i++) guess[i] -= lastderiv[i+size];

            for(i=0,k=0;i<size;i++) {
                for(j=0;j<i;j++){
                    curvature[k++] += sums * guess[i] * guess[j];
                }
                curvature[k++] += sums * guess[i] * guess[i];
            }
            //        if ((shrinkcount & 7) == 7) {for(i=0;i<size;i++) lastderiv[i] = -lastderiv[i];}

            memcpy(guess,lastderiv+size,sizeof(double)*size);
            //   printf("curvature:");

            modif = 0.0f;

            for(i=0,k=0;i<size;i++) {
                guess[i] += curvature[k] * lastderiv[i];//printf("%e\t",curvature[k] );
                for(j=i+1,k++;j<size;j++,k++){
                    guess[j] += curvature[k] * lastderiv[i];//printf("%e\t",curvature[k] );
                    guess[i] += curvature[k] * lastderiv[j];

                }
                modif += (fabs(guess[i]) > 1.2e-100) ? fabs((guess[i] - lastderiv[i+size])/ guess[i]) : 0.0f ;
            }//printf("\n");


            return( modif );
        }else if ((sump < 0.9f * sums)&&(alpha &0x7FFFFFFF)) { // step too short
            // scaling inverse hessian by 4
            // sums is the denominator!
            //     printf("Expand!\n");
            alpha--;

            sums = 3.0f / fabs(sums);
            shrinkcount =0;
            for(i=0;i<size;i++) lastderiv[i+size] -= guess[i];

            for(i=0,k=0;i<size;i++) {
                curvature[k] += sums * (lastderiv[i+size] * lastderiv[i+size]);// printf("%f\t", curvature[k]);
                for(j=i+1,k++;j<size;j++,k++){
                    curvature[k] += sums * (lastderiv[j+size] * lastderiv[i+size]); //printf("%f\t", curvature[k]);
                }
            }


            memcpy(lastderiv + size ,guess,sizeof(double)*size);
            //    printf("curvature:");
            for(i=0,k=0;i<size;i++) {
                guess[i] += curvature[k] * deriv[i];//printf("%e\t",curvature[k] );
                for(j=i+1,k++;j<size;j++,k++){
                    guess[j] += curvature[k] * deriv[i];
                    guess[i] += curvature[k] * deriv[j];//printf("%e\t",curvature[k] );
                }

            }// printf("\n");
            lastvalue = value;
            memcpy(lastderiv,deriv,sizeof(double)*size);

            return(ExCo<double>::mkMaximum());
        }else { shrinkcount =0;
            // if condition1 fails, the step is too long: it is impossible to update the hessian
            // if condition2 fails, the step is too short: it is impossible to update the hessian
         //   printf("AlRight!\n");

            for(i=0;i<size;i++) {
                lastderiv[i+size] -= guess[i]; // diff
                lastderiv[i+3*size] = deriv[i] - lastderiv[i];
            }

            // lastderiv[size*2+x] = B^-1 (deriv-lastderiv)
            k=0;
            lastderiv[size*2] = curv.data[k++] * lastderiv[3*size];
            for(i=1;i<size;i++) {
                lastderiv[size*2] += curv.data[k] * lastderiv[i+3*size];
                lastderiv[i+size*2] = curv.data[k++] * lastderiv[3*size];
                for(j=1;j<i;j++){
                    lastderiv[j+size*2] += curv.data[k] * lastderiv[i+3*size];
                    lastderiv[i+size*2] += curv.data[k++] * lastderiv[j+3*size];
                }
                lastderiv[i+size*2] += curv.data[k++] * lastderiv[i+3*size];
            }


            sums -= sump;

            sump = -sums;
            for(i=0;i<size;i++) sump += lastderiv[i+size*2] * lastderiv[i+3*size];
            sums = -1.0f / sums;
            sump *= sums;
        //     printf("K2 = %f\tK1 = %f\n", sump,sums);
            for(i=0,k=0;i<size;i++) {
                for(j=0;j<i;j++){
                    curv.data[k++] += sums *(sump * lastderiv[i+size]* lastderiv[j+size] - (lastderiv[i+size]* lastderiv[j+size*2] + lastderiv[i+size*2]* lastderiv[j+size]));
                }
                curv.data[k++] += sums * (lastderiv[i+size] * (sump * lastderiv[i+size] - 2.0f * lastderiv[i+size*2]));
            }

        }
    }

    // has updated curvature by now

    // renormalize hessian to match provided double derivatives


    //
    modif = (fabs(value) > 1.2e-100) ? fabs((value - lastvalue) / (value)) : 0.0f;
    memcpy(lastderiv + size ,guess,sizeof(double)*size);
//   printf("curvature:");

    double fact = pow(0.5f, (double)(alpha&0x7FFFFFFF));
    if (alpha & 0x80000000) fact = -fact;

    expdiff = 0.0f;
    for(i=0,k=0;i<size;i++) {
        for(j=0;j<i;j++){
            guess[j] += curv.data[k] * fact * deriv[i];
            guess[i] += curv.data[k] * fact * deriv[j];//printf("%e\t",curvature[k] );
            expdiff += curv.data[k++] * deriv[i] * deriv[j];
        }
        guess[i] += curv.data[k] * fact * deriv[i];//printf("%e\t",curvature[k] );
        expdiff += 0.5f * curv.data[k++] * deriv[i] * deriv[i];
        //modif += (fabs(guess[i]) > 1.2e-100) ? fabs((guess[i] - lastderiv[i+size])/ guess[i]) : 0.0f ;
    } //printf("\n");
    expdiff *= 2.0f;

    lastvalue = value;
    memcpy(lastderiv,deriv,sizeof(double)*size);
//    printf("%e\n",modif /size);

    return(modif);
}*/

	double CurvatureSearchScope::updateDescent(const double &value, double * guess , double const * const deriv){ // returns a safe bound on the parameter changes;
		unsigned int i,j,k;
		double sums,  sump, modif;
		if (lastderiv == NULL) {
            lastderiv = new double[size*4];
  		    if (!ExOp::isValid(value)) {fprintf(stderr, "initial guess has function value NaN, provide a better guess...\n"); exit(1);}
            modif = ExCo<double>::mkMaximum();
        }else if (!ExOp::isValid(value)){ // not-a-number eh?
            memcpy(guess,lastderiv+size,sizeof(double)*size);
            modif = 0.0f;
            for(i=0,k=0;i<size;i++) {
                curvature[k] *= 0.0625;
                guess[i] -= curvature[k] * lastderiv[i];
                for(j=i+1,k++;j<size;j++,k++){
                    curvature[k] *= 0.0625;
                    guess[j] -= curvature[k] * lastderiv[i];
                    guess[i] -= curvature[k] * lastderiv[j];

                }
                modif += (fabs(guess[i]) > 1.2e-100) ? fabs((guess[i] - lastderiv[i+size])/ guess[i]) : 0.0f ;
            }
            return modif;
        }else{

			// check Wolfe Conditions!

			// value - lastvalue = 0.0001(guess - oldguess) \cdot oldderiv
			// (guess - oldguess) \cdot deriv = 0.9 (guess - oldguess) \cdot oldderiv
            sums =0; sump =0;
            for(i=0;i<size;i++) {sums += lastderiv[i] * (guess[i] - lastderiv[i+size]);
                                     sump += deriv[i] * (guess[i] - lastderiv[i+size]);}
        //    printf("wolfe 1: %e <= %e  %c\n",  value - lastvalue, 0.0001f * sums, (value - lastvalue <= 0.0001f * sums) ? 'Y' : 'N' );
       //     printf("wolfe 2: %e >= %e  %c\n", sump, 0.9f * sums , (sump>= 0.9f * sums) ? 'Y' : 'N' );
       //      printf("wolfe 1: %e <= %e  %c\n", value - lastvalue, 0.0001f * sums, (value - lastvalue <= 0.0001f * sums) ? 'Y' : 'N' );
       //     printf("wolfe 2: %e >= %e  %c\n", sump, 0.9f * sums , (sump>= 0.9f * sums) ? 'Y' : 'N' );
        if (!ExOp::isValid(value)){
            memcpy(guess,lastderiv+size,sizeof(double)*size);
            modif = 0.0f;
            for(i=0,k=0;i<size;i++) {
                curvature[k] *= 0.0625f;
                guess[i] += curvature[k] * lastderiv[i];//printf("%e\t",curvature[k] );
                for(j=i+1,k++;j<size;j++,k++){
                    curvature[k] *= 0.0625f;
                    guess[j] += curvature[k] * lastderiv[i];//printf("%e\t",curvature[k] );
                    guess[i] += curvature[k] * lastderiv[j];

                }
                modif += (fabs(guess[i]) > 1.2e-100) ? fabs((guess[i] - lastderiv[i+size])/ guess[i]) : 0.0f ;
            }
            return modif;
        }
        if (value - lastvalue > 0.0001f * sums){  // step too long
            // scaling inverse hessian by 0.25
			// sums is the denominator!

			//   printf("desc Shrink!\n");
            sums = -0.75f / fabs(sums);
            sump=0;
            for(i=0;i<size;i++) guess[i] -= lastderiv[i+size];

            for(i=0,k=0;i<size;i++) {
                curvature[k] += sums * guess[i] * guess[i];// printf("%f\t", curvature[k]);
                for(j=i+1,k++;j<size;j++,k++){
                    curvature[k] += sums * guess[i] * guess[j]; //printf("%f\t", curvature[k]);
                }
            }

            memcpy(guess,lastderiv+size,sizeof(double)*size);
			//     printf("curvature:");
            for(i=0,k=0;i<size;i++) {
                guess[i] -= curvature[k] * lastderiv[i];//printf("%e\t",curvature[k] );
                for(j=i+1,k++;j<size;j++,k++){
                    guess[j] -= curvature[k] * lastderiv[i];//printf("%e\t",curvature[k] );
                    guess[i] -= curvature[k] * lastderiv[j];
                }
            }//printf("\n");


            return(ExCo<double>::mkMaximum());
        }else if (sump < 0.9f * sums){ // step too short
            // scaling inverse hessian by 4
            // sums is the denominator!
			//     printf("desc Expand!\n");
            sums = 3.0f / fabs(sums);

            for(i=0;i<size;i++) lastderiv[i+size] -= guess[i];

            for(i=0,k=0;i<size;i++) {
                curvature[k] += sums * (lastderiv[i+size] * lastderiv[i+size]);// printf("%f\t", curvature[k]);
                for(j=i+1,k++;j<size;j++,k++){
                    curvature[k] += sums * (lastderiv[j+size] * lastderiv[i+size]); //printf("%f\t", curvature[k]);
                }
            }


			memcpy(lastderiv + size ,guess,sizeof(double)*size);
			//    printf("curvature:");
			for(i=0,k=0;i<size;i++) {
				guess[i] -= curvature[k] * deriv[i];//printf("%e\t",curvature[k] );
				for(j=i+1,k++;j<size;j++,k++){
					guess[j] -= curvature[k] * deriv[i];
					guess[i] -= curvature[k] * deriv[j];//printf("%e\t",curvature[k] );
				}

			} //printf("\n");
			lastvalue = value;
			memcpy(lastderiv,deriv,sizeof(double)*size);

			return(ExCo<double>::mkMaximum());
        }else {
            // if condition1 fails, the step is too long: it is impossible to update the hessian
			// if condition2 fails, the step is too short: it is impossible to update the hessian

//printf("desc Right!\n");
        for(i=0;i<size;i++) lastderiv[i+size] -= guess[i];
        i=0;k=0;

        lastderiv[i+size*2] = curvature[k] * deriv[i];

        for(j=i+1,k++;j<size;j++,k++){
            lastderiv[j+size*2] = curvature[k] * deriv[i];
            lastderiv[i+size*2] += curvature[k] * deriv[j];
   //         printf("%f\t", curvature[k]);
        }
        for(i++;i<size;i++) {
            lastderiv[i+size*2] += curvature[k] * deriv[i];
            for(j=i+1,k++;j<size;j++,k++){
                lastderiv[j+size*2] += curvature[k] * deriv[i];
                lastderiv[i+size*2] += curvature[k] * deriv[j];
            }
        }

        sums -= sump;
        sump = lastderiv[size*2] * (deriv[0] - lastderiv[0]);
        for(i=1;i<size;i++) sump += lastderiv[i+size*2] * (deriv[i] - lastderiv[i]);
        sums = 1.0f / sums;
        sump *= sums;
   //     printf("K2 = %f\tK1 = %f\n", sump,sums);
        for(i=0,k=0;i<size;i++) {
            curvature[k] += sums * (lastderiv[i+size] * (sump * lastderiv[i+size] - 2.0f * lastderiv[i+size*2]));
            for(j=i+1,k++;j<size;j++,k++){
                curvature[k] += sums *(sump * lastderiv[i+size]* lastderiv[j+size] - (lastderiv[i+size]* lastderiv[j+size*2] + lastderiv[i+size*2]* lastderiv[j+size]));
            }
        }

		}
            modif = (fabs(value) > 1.2e-100) ? fabs((value - lastvalue) / (value)) : 0.0f;
        }

        memcpy(lastderiv + size ,guess,sizeof(double)*size);
     //    printf("curvature:");
        for(i=0,k=0;i<size;i++) {
            guess[i] -= curvature[k] * deriv[i];//printf("%e\t",curvature[k] );
            for(j=i+1,k++;j<size;j++,k++){
                guess[j] -= curvature[k] * deriv[i];
                guess[i] -= curvature[k] * deriv[j];//printf("%e\t",curvature[k] );
            }
         //   modif += (fabs(guess[i]) > 1.2e-100) ? fabs((guess[i] - lastderiv[i+size])/ guess[i]) : 0.0f ;
        } //printf("\n");
        lastvalue = value;
        memcpy(lastderiv,deriv,sizeof(double)*size);
    //    printf("%e\n",modif /size);

		return(modif);
	}

void CurvatureSearchScope::wrFinalGuess(double* guess)const{if (lastderiv != NULL) memcpy(guess,lastderiv + size,sizeof(double)*size);}
double CurvatureSearchScope::checkDerivative(const double &value, double * guess , double const * const deriv){
    if (lastderiv == NULL){ // first time, dont update alpha
        lastderiv = new double[size+2];
        guess[0] += curvature[0];
        lastderiv[size] =0;
    }else{
        int which = (int)lastderiv[size];
        printf("%i: %f%% error, %e and %e predicted, %e computed.\n", which, ((lastderiv[which] + deriv[which]) * curvature[0] * 0.5f / (value - lastvalue)) -1.0f, lastderiv[which], deriv[which],  (value - lastvalue) / curvature[0]);
        which = (which+1) % size;
        guess[which] += curvature[0];
        lastderiv[size] = which;
    }
    lastvalue = value;
    memcpy(lastderiv,deriv,sizeof(double)*size);
    return(curvature[0]);
}

double FuzzyLineSearch::initGuesses(double cur_best, double newguess){
    ExOp::toZero(value);
    guess[0] = cur_best;
    guess[1] = (cur_best + newguess) * 0.5;
    guess[2] = newguess;
}
double FuzzyLineSearch::solveArgmax(double& cur_best, double &nextguess) const{ // find the best of the 3 guess, and return next guess
    double discrim = (value[0] + value[2] - value[1] * 2.0) * 4.0;
    uint32_t bindex;
    if (!ExOp::isValid(discrim)){ // GOT NA values! in Therory value[0] should not be NA
        if ((ExOp::isValid(value[1]))&&(value[1] > value[0])){
            bindex =1;
            nextguess = guess[1] * 1.5 - guess[2] * 0.5;
        }else if ((ExOp::isValid(value[2]))&&(value[2] > value[0])){
            bindex =2;
            nextguess = guess[2] * 1.5 - guess[0] * 0.5;
        }else{
            bindex =0;
            nextguess = guess[0] *  0.75 + guess[1] * 0.25;
        }
    }else{
        double fout;

        double offset;
        if ((discrim >= 0.0)||(!ExOp::isValid(offset = (value[2] - value[0]) / discrim))){ // incredible, wrong curvature or overflow
            if  (value[0] > value[2]){
                bindex=0;
                nextguess = guess[0] * 16.5 - guess[2] * 15.5;
            }else{
                bindex=2;
                nextguess = guess[2] * 16.5 - guess[0] * 15.5;
            }
        }else{ // find maximum offset
            if (offset > 0.0){
                if (offset > 16.0) offset = 16.0;
                bindex = (value[0] > value[1]) ? 0 : 1;
            }else{
                if (offset < -16.0) offset = -16.0;
                bindex = (value[2] > value[1]) ? 2 : 1;
            }
            nextguess = guess[1] - (guess[2] - guess[0]) * offset;
        }
    //    printf("discr %e %e ->%e\n",discrim,offset, nextguess);
    //    printf("%e %e %e -> %i  %e\n",value[0],value[1],value[2],bindex, nextguess);
    }
    cur_best = guess[bindex];
return value[bindex];}

UnknownScope::UnknownScope(double freq){scp[1] = freq; scp[2] = 0.0f; oldLL = -DBL_MAX; memset(scp+3,'\0',sizeof(double)*5);}

// p out = 1 / (1 + exp(D - x))

bool UnknownScope::EMprocess(double* LLs, double &total_LL, uint32_t length, Tuple<double,2u> *outlier_prior){
    total_LL =0.0;
    double tmp;
    uint32_t i;
    if (scp[1] == 0.0){
        if (outlier_prior == NULL){
            for(i=0;i<length;i++){
                total_LL += LLs[i];
                LLs[i] = 1.0f;
            }
        }else{
            for(i=0;i<length;i++){
                tmp = LLs[i];
                LLs[i] = (1.0 + outlier_prior[i][0]) / (1.0 + outlier_prior[i][0] + outlier_prior[i][1]);
                total_LL += LLs[i] * tmp;
            }
        }
        if (oldLL < total_LL){
            oldLL = total_LL;
            return true;
        }else return false;
    }

    double pivots[6];
    double stat[6];

    memset(scp+3,'\0',sizeof(double)*5);

    double logittarget = log(scp[1]) - log(1.0 - scp[1]);
    i=0;


    pivots[1] = pivots[i] = LLs[i];
    stat[0] = LLs[i];
    stat[1] = LLs[i] * LLs[i];
    for(i++;i<length;i++){
        stat[0] += LLs[i];
        stat[1] += LLs[i] * LLs[i];
        if (pivots[1] < LLs[i]) pivots[1] =LLs[i];
        else if (pivots[0] > LLs[i]) pivots[0] =LLs[i];
    }

    double priorbase = 0.0;
    if (outlier_prior != NULL){
        for(i=0;i<length;i++){
            priorbase  += outlier_prior[i][0] / (1.0 + outlier_prior[i][0] + outlier_prior[i][1]);
        }
    }

    //printf("ss %e,%e\n", stat[0],stat[1]);
    stat[0] /= length;
    stat[1] = sqrt((stat[1] / length) - stat[0] * stat[0]);
    //printf("%e +- %e\n", stat[0],stat[1]);
    //printf("logit target %e\n", logittarget);

    //printf("F[%e] = %e\n", scp[1], sliwa_gsl_cdf_ugaussian_Pinv(scp[1]));
    if (ExOp::isValid(scp[5])) pivots[2] = stat[0] + stat[1] * sliwa_gsl_cdf_ugaussian_Pinv(scp[1]) ; // aim for half of all.. "mean" LL    fr = 1 / (1 + exp(LL - x) ) ,  logit(fr) = LL - x, x = LL - logit(fr)
    else {pivots[2] = stat[0] * scp[1];}
    memset(stat,'\0',sizeof(double)*3);
    if (outlier_prior != NULL){
        for(i=0;i<length;i++){
            tmp = (1.0 + outlier_prior[i][0] + outlier_prior[i][1]);
            stat[0] += logistic(pivots[2] - LLs[i]) / tmp;
            stat[1] += logistic(pivots[0] - LLs[i]) / tmp;
            stat[2] += logistic(pivots[1] - LLs[i]) / tmp;
        }
    }else{
        for(i=0;i<length;i++){
            stat[0] += logistic(pivots[2] - LLs[i]);
            stat[1] += logistic(pivots[0] - LLs[i]);
            stat[2] += logistic(pivots[1] - LLs[i]);
        }
    }
    pivots[5] = (priorbase + stat[0]) / length; pivots[3] = (priorbase +stat[1]) / length; pivots[4] = (priorbase + stat[2]) / length;
   // printf("Init search F[%e %e %e] = %e %e %e\n",pivots[0], pivots[2], pivots[1], pivots[3], stat[0],pivots[4]);
    pivots[3] -= scp[1];
    pivots[4] -= scp[1];
    pivots[5] -= scp[1];
    /*for(uint32_t k=0;k<100;k++){ // ez secant method
        if (pivots[5] < 0){ pivots[0] = pivots[2]; pivots[3] = pivots[5];}
        else {pivots[1] = pivots[2]; pivots[4] = pivots[5];}
        scp[0] = (pivots[1] - pivots[0]) * pivots[3] / (pivots[3] - pivots[4]) + pivots[0];
        stat[0] = 0.0;
        for(i=0;i<length;i++) stat[0] += logistic(scp[0] - LLs[i]);
        stat[0] /= length;
        printf("F[%e] = %e\n",scp[0], stat[0]);
        pivots[2] = scp[0];
        pivots[5] = stat[0] - scp[1];
    }*/

    for(uint32_t k=0;k<100;k++){ // quadratic
        stat[0] = pivots[3] / ((pivots[0] - pivots[1]) * (pivots[0] - pivots[2]));
        stat[1] = pivots[4] / ((pivots[1] - pivots[0]) * (pivots[1] - pivots[2]));
        stat[2] = pivots[5] / ((pivots[2] - pivots[0]) * (pivots[2] - pivots[1]));
        stat[3] = stat[0] + stat[1] + stat[2];
        stat[4] = stat[0] *(pivots[1] + pivots[2])+ stat[1] *(pivots[0] + pivots[2])+ stat[2] *(pivots[0] + pivots[1]);
        stat[4] /= stat[3] * 2.0;
        stat[5] = stat[0] * pivots[1] * pivots[2] + stat[1] * pivots[0] * pivots[2] + stat[2] * pivots[0] * pivots[1];
        stat[0] = sqrt(stat[4] * stat[4] - stat[5] / stat[3]);

        scp[0] = (stat[3] < 0.0) ? stat[4]  - stat[0] : stat[4] + stat[0];
        if (!ExOp::isValid(scp[0])){ // has gone mad, use secant method instead
            if (pivots[5] < 0.0f){ pivots[0] = pivots[2]; pivots[3] = pivots[5];}
            else {pivots[1] = pivots[2]; pivots[4] = pivots[5];}
            scp[0] = (pivots[1] - pivots[0]) * pivots[3] / (pivots[3] - pivots[4]) + pivots[0];
            if (!ExOp::isValid(scp[0])) {scp[0] = pivots[1]; break;} // probably converged already (or lost)
        }else{
            if (scp[0] < pivots[2]) { pivots[1] = pivots[2]; pivots[4] = pivots[5];}
            else {pivots[0] = pivots[2]; pivots[3] = pivots[5];}
        }
        stat[0] = 0;

        if (outlier_prior != NULL) for(i=0;i<length;i++) stat[0] += logistic(scp[0] - LLs[i]) / (1.0 + outlier_prior[i][0] + outlier_prior[i][1]);
        else for(i=0;i<length;i++) stat[0] += logistic(scp[0] - LLs[i]);


        stat[0] = (priorbase + stat[0]) / length;
      //  printf("F[%e] = %e\n",scp[0], stat[0] - scp[1]);
        pivots[2] = scp[0];
        pivots[5] = stat[0] - scp[1];
    }
    if (outlier_prior != NULL){
        for(i=0;i<length;i++){
            scp[4] = LLs[i];
            LLs[i] = (logistic(LLs[i] - scp[0]) + outlier_prior[i][0]) / (1.0 + outlier_prior[i][0] + outlier_prior[i][1]);
            total_LL += LLs[i] * scp[4];
        }
    }else{
        for(i=0;i<length;i++){
            scp[4] = LLs[i];
            LLs[i] = logistic(LLs[i] - scp[0]);
            total_LL += LLs[i] * scp[4];
        }
    }

    if (oldLL < total_LL){
        oldLL = total_LL;
        return true;
    }else return false;
}

double UnknownScope::EMregist(double LL){
    int i;
    double sum,tmp,tmp2;
    double pix[8];

    if (scp[2] == 0.0f){ // passive initial guess
        if (ExOp::isValid(LL)) {
            scp[3] += 1.0f;
            scp[4] += LL;
            scp[5] += LL * LL;
            }

        return 1.0;
    }else { // close approach new strat F(x) dF(x) ddF(x)  F(x+step)
        sum = exp(LL - scp[0]);
        tmp = 1.0f + sum;
        pix[0] = logistic(scp[0] - LL);
        tmp2 = 2.0f + sum + (1.0f / sum);
        pix[1] = 1.0f / tmp2;
        pix[2] = 1.0f / (tmp2 * tmp);
        pix[3] = 1.0f / (tmp2 * tmp * tmp);
        if ((ExCo<double>::isValid(pix[0]))&&(ExCo<double>::isValid(pix[1]))&&(ExCo<double>::isValid(pix[2]))&&(ExCo<double>::isValid(pix[3]))){
            scp[3] += 1.0f;
            scp[4] += pix[0]; // real
            scp[5] += pix[1];
            scp[6] += pix[2];
            scp[7] += pix[3];
        }
        return pix[0];
    }

}

bool UnknownScope::EMfinit(double &total_LL){
    double coef[4];
    if (scp[2] > 0.125f){ // far approach
        //		printf("Far Step:\n");
        //		printf("%e\t%e\t%e\t%e\n", (scp[4]-scp[5])/(scp[3]), (scp[6]-scp[7])/(scp[3]),(scp[6]+scp[7])/(scp[3]),(scp[4]+scp[5])/(scp[3]));

        coef[0] =  (-scp[4] + 9.0f * scp[6]) - 8.0f * scp[3] * scp[1];
        coef[1] =  ((scp[5]/-3.0f) + 9.0f * scp[7]) ;
        coef[2] =  (scp[4] - scp[6]);
        coef[3] =  ((scp[5]/3.0f) - scp[7]);
        //		printf("%e\t%e\t%e\t%e\n", coef[0], coef[1], coef[2], coef[3]);
        double shift = CubicRealRoot(coef, false);
        //printf("Far Step: F[%e] = %e!\n",exp(scp[0]), scp[6]/(scp[3]));
        //		printf("%e\t%e\n", shift, scp[2]);

        shift *=antioss;
        if (!ExCo<double>::isValid(shift)){
            scp[0] += (scp[6] > scp[3]* scp[1] ? -3.0f : 3.0f)* scp[2];

        }else if (fabs(shift) > 3.0f){ // outside!
            scp[0] += (shift < 0 ? -3.0f : 3.0f)* scp[2];
            scp[2] *= 4.0f;

        }else{
            scp[0] += shift * scp[2];
            scp[2] *= fabs(scp[1] -(scp[5]+ scp[6])/(scp[3]));
        }
        //		printf("end: %e\t    %e\n", exp(scp[0]), scp[2]);
        relerr = log(scp[6]) - log(scp[3] * scp[1]);

        if (relerr < 1.0f){ // dont update if too much data is in unknown class
            total_LL /= (1.0f - scp[1] * exp(relerr));
        } else {scp[0] -= 3.0;antioss *=0.9f; }
    }else if (scp[2] == 0.0f){
        printf("ss %e,%e,%e\n", scp[3],scp[4],scp[5]);
        scp[4] /= scp[3];
        scp[5] = sqrt((scp[5] / scp[3]) - scp[0] * scp[0]);
        printf("%e +- %e\n", scp[4],scp[5]);

        printf("F[%e] = %e\n", scp[1], sliwa_gsl_cdf_ugaussian_Pinv(scp[1]));
        if (ExOp::isValid(scp[5])) scp[0] = scp[4] - scp[5] * sliwa_gsl_cdf_ugaussian_Pinv(scp[1]) ; // aim for half of all.. "mean" LL    fr = 1 / (1 + exp(LL - x) ) ,  logit(fr) = LL - x, x = LL - logit(fr)
        else {scp[0] = scp[4] * scp[1];}
        printf("new outlierLL %e\n", scp[0]);
        scp[2] = 20.0f;
        relerr = 20.0f;
        antioss = 1.0f;
        total_LL = scp[4];
        oldLL = scp[4]; // likelihood with all data, should be lower than partial data
        memset(scp+3,'\0',sizeof(double)*5);
        return true;
    }else{ // close approach
        double tmp;
        //	printf("Close Step:\n");
        relerr = log(scp[4]) - log(scp[1] * scp[3]);
        if ((ExCo<double>::isValid(scp[6]))&&(ExCo<double>::isValid(scp[7]))&&(scp[6] != 0.0f)&&(scp[7] != 0.0f)&&(fabs(1.0f / scp[6]) != 0.0f) &&(fabs(1.0f / scp[7]) != 0.0f) ){

            coef[0] = scp[4] - scp[1] * scp[3];
            coef[1] = scp[5];
            coef[2] = scp[6];
            coef[3] = scp[7];
            //		poly[2] = 2 * statbuf[3] + statbuf[2];
            //		poly[3] = 6 * statbuf[4] + 6 * statbuf[3] + statbuf[2];
            //		printf("%e\t%e\t%e\t%e\n",coef[0], coef[1], coef[2] , coef[3]  );

            tmp = -CubicRealRoot(coef,false);

            //		printf("%e poly eval\n", coef[0] -tmp *(coef[1] -tmp * (coef[2] -tmp * coef[3])));

            if ((!ExCo<double>::isValid(tmp))||(fabs(tmp) > 2.0f)){
                tmp = (scp[4] - scp[1] * scp[3]) / scp[5];
                //				printf("%e\t%e\n",(scp[4] - scp[1] * scp[3]) ,  scp[5]);
                if (!ExCo<double>::isValid(tmp)) tmp = 2.0f;
            }

        }else if (scp[3] > 0.0f){

            tmp = (scp[4] - scp[1] * scp[3]) / scp[5];
            //		printf("wrongwrong F[%e] = %e ! %e\n",tmp ,scp[3],scp[4]);
            //		printf("%e\t%e\n",(scp[4] - scp[1] * scp[3]) ,  scp[5]);
            if (!ExCo<double>::isValid(tmp)) tmp = 2.0f;


        } else {tmp = 2.0f; scp[0] -= 3.0;antioss *=0.9f; relerr = 4.0f;}
        //printf("Near Step: F[%e] = %e !\n",exp(scp[0]), scp[4] / scp[3]);

        tmp *=antioss;
        if (fabs(tmp) >= 2.0f){
            scp[0] +=  (tmp < 0.0f) ? 2.0 : -2.0;
            scp[2] *= 1.25f;
            if (scp[2] > 20.0f) scp[2] = 0.0f;
        }else{
            scp[2] *= fabs(tmp) /2.0f;
            scp[0] -= tmp;
        }
        //		printf("step = %e, new guess %e\n", tmp,exp(scp[0]));

        if ((relerr < 1.0f)&&(ExCo<double>::isValid(relerr))){ // dont update if too much data is in unknown class
            total_LL /= (1.0f - scp[1] * exp(relerr));
        }else {scp[0] -= 3.0;antioss *=0.9f;} // want to converge from below
    }
    memset(scp+3,'\0',sizeof(double)*5);
    printf("new outlierLL %e\n", scp[0]);
return true;}

// F(x) = F(x_0) + (x-x_0) * dF(x_0) +- K * (x-x_0)^2  -> F(x) \sim N (F(x_0) + (x-x_0) * dF(x_0), K*(x-x_0)^2)

// x' = x_0 -



Dictionary::Dictionary(): flags(0){}
Dictionary::~Dictionary(){this->toMemfree();}
Dictionary& Dictionary::toMemfree(){
    for(uint32_t i=0;i< entries.getSize();i++) delete[](entries[i]);
    entries.toMemfree();
    word_links.toMemfree();
    hashes.toMemfree();
return *this;}
Dictionary& Dictionary::toMemmove(Dictionary& source){
    entries.toMemmove(source.entries);
    word_links.toMemmove(source.word_links);
    hashes.toMemmove(source.hashes);
return *this;}

uint32_t Dictionary::find(const char* what) const{
    return findWord(what, strlen(what));
    /*uint32_t key;
    uint32_t selhash = 0;
    if (strlen(what) <3){
        if (strlen(what) == 2) ((((uint32_t)what[0]) << 8) | what[1]) << 16;
        else key = (((uint32_t)what[0]) << 24);
    }else key = (((((((uint32_t)what[0]) << 8) | what[1]) << 8) | what[2])) | what[3];
    ite = hashes[selhash].find(key);
    while(ite != 0xFFFFFFFF){
        selhash =  hashes[selhash].deref(ite);

    }
    return 0xFFFFFFFF;*/
}

uint32_t Dictionary::findWord(const char* word, int word_lenght) const{
	//printf("find lenght: %i\n", word_lenght);

	int proc = 0;
	uint32_t key;
	uint32_t selhash = 0;
	uint32_t ite;
    if (hashes.getSize() == 0) return 0xFFFFFFFF;
	while(word_lenght - proc >= 4){
		key = (((int)word[proc]) << 24) | (((int)word[proc | 1]) << 16) | (((int)word[proc | 2]) << 8) | ((int)word[proc | 3]);
		ite = hashes[selhash].find(key);
		if (ite == 0xFFFFFFFF) return 0xFFFFFFFF;
		selhash = hashes[selhash].deref(ite);
		//printf("%i %i tmp\n", selhash, proc);
		//hashes[selhash].show();
		proc+= 4;
	}
	if (word_lenght - proc >= 2){
		key = (((int)word[proc | 1]) << 16) | ((word_lenght - proc == 3) ?  (((int)word[proc | 2]) << 8) : 0);
	} else key = 0;
	if (word_lenght - proc > 0) key |= (((int)word[proc]) << 24);
	ite = hashes[selhash].find(key);
	//if (ite == 0xFFFFFFFF) printf("cant find %s\n", word);
	//else printf("found %s at %X\n", word,hashes[selhash].deref(ite));
	return (ite == 0xFFFFFFFF) ? 0xFFFFFFFF : hashes[selhash].deref(ite);
}
uint32_t Dictionary::insertWord(const char* word, int word_lenght){
	//printf("inserting word %s\n", word); fflush(stdout);
	uint32_t fout = word_links.getSize();
	word_links.push_back();
//	word_links[fout].k = new char[word_lenght+1];
//	memcpy(word_links[fout].k,word, word_lenght);
//	word_links[fout].k[word_lenght] = '\0';
	int proc = 0;
	unsigned int ite;
	uint32_t selhash = 0;
	if (hashes.getSize() ==0 ) hashes.setSize(1);
	uint32_t key;
	while(word_lenght - proc >= 4){
		key = (((int)word[proc]) << 24) | (((int)word[proc | 1]) << 16) | (((int)word[proc | 2]) << 8) | ((int)word[proc | 3]);
		if ((ite = hashes[selhash].find(key)) != 0xFFFFFFFF) {
            selhash =  hashes[selhash].deref(ite);
          //  printf("proc%i got to %i\n", proc, selhash);
		}else{
		   // printf("creating new hash for %X\n", key);
			hashes[selhash][key] = hashes.getSize();
			hashes[selhash].keysort();
			selhash = hashes.getSize();
			hashes.push_back();
		}
		proc+= 4;
	}

	if (word_lenght - proc >= 2){
		key = (((int)word[proc | 1]) << 16) | ((word_lenght - proc == 3) ?  (((int)word[proc | 2]) << 8) : 0);
	} else key = 0;
	if (word_lenght - proc > 0) key |= (((int)word[proc]) << 24);
	ite = hashes[selhash].find(key);
	if (ite != 0xFFFFFFFF) {fprintf(stderr, "DICO error: trying to insert word that already exist!\n"); exit(1);}
//	printf("final entry within %X, (index %i)", key, selhash);
	hashes[selhash][key] = fout;
	hashes[selhash].keysort();
	return fout;
}
Tuple<uint32_t, 4u> Dictionary::findPrefix(const char* prefix, int word_lenght) const{
	Tuple<uint32_t, 4u> fout;
	int proc = 0;
	uint32_t key;
	uint32_t selhash = 0;
	uint32_t ite;
	if (hashes.getSize() ==0 ) {fout[0]= 0xFFFFFFFF; return fout;}

	while(word_lenght - proc >= 4){
		key = (((int)prefix[proc]) << 24) | (((int)prefix[proc | 1]) << 16) | (((int)prefix[proc | 2]) << 8) | ((int)prefix[proc | 3]);
		ite = hashes[selhash].find(key);
		if (ite == 0xFFFFFFFF) {fout[0]= 0xFFFFFFFF; return fout;}
		selhash = hashes[selhash].deref(ite);
		proc+= 4;
	}
	if (word_lenght - proc >= 2){
		key = (((int)prefix[proc | 1]) << 16) | ((word_lenght - proc == 3) ?  (((int)prefix[proc | 2]) << 8) : 0);
	} else key = 0;
	if (word_lenght - proc > 0) key |= (((int)prefix[proc]) << 24);

	fout[0] = 0;
	fout[1] = hashes[selhash].getSize()-1;
	//printf("binary search!\n");

	while(fout[0] != fout[1]){
		if (hashes[selhash].deref_key((fout[0]+fout[1]) >> 1) >= key) fout[1] = (fout[0]+fout[1]) >> 1;
		else fout[0] = ((fout[0]+fout[1]) >> 1) +1;
	}


	fout[2] = key;
	fout[3] = (word_lenght == proc) ? 0 : 0xFFFFFFFF << ((4 - word_lenght + proc) * 8);

	//printf("%X and %X xor : %X & %X\n", key, hashes[selhash].deref_key(fout[1]), (hashes[selhash].deref_key(fout[1]) ^ key), fout[3]);

	ite = hashes[selhash].deref_key(fout[1]) ^ key;



	fout[0] = (ite & fout[3]) ? 0xFFFFFFFF : selhash;

	//if (ite == 0xFFFFFFFF) printf("cant find %s\n", word);
	//else printf("found %s at %X\n", word,hashes[selhash].deref(ite));
	return fout;

}
void Dictionary::addEntry(const char* str, bool allow_duplicate){
	uint32_t i = strlen(str);
	char* newstr = new char[i+1];
	memcpy(newstr, str, i+1);
	uint32_t word_index = this->findWord(newstr, i);
	if (word_index != 0xFFFFFFFF) {
            if (!allow_duplicate) {fprintf(stderr, "warning: dictionary already contains %s\n", str); return;}
    }
	memcpy(newstr, str, i+1);
	entries.push_back(newstr);

	int c;
	if (flags & 1){
		for(c=0;newstr[c] != '\0';c++){
			if (((newstr[c] & 0xC0) == 0x40)&&(newstr[c] <= 'Z')) newstr[c] |= 0x20;
		}
	}


	const char* c_start = newstr;
	const char* c_end = newstr;

	while(true){
		while(true){
			if ((*c_end == '\0')||(*c_end == ' ')||(*c_end == '\t')) break;
			c_end++;
		}
		if (c_end != c_start){
			word_index = this->findWord(c_start, (int)(c_end - c_start));
			if (word_index == 0xFFFFFFFF) word_index = insertWord(c_start, (int)(c_end - c_start));
			for(i=0;i<word_links[word_index].getSize();i++) if (word_links[word_index][i] == entries.getSize()-1) break;
			if (i==word_links[word_index].getSize()) word_links[word_index].push_back(entries.getSize()-1);
		}
		if (*c_end == '\0') break;
		c_end++;
		c_start = c_end;
	}

	if (flags & 1) memcpy(newstr, str, i+1);
}
void Dictionary::registerSubstrings(int maxlen){
	unsigned int i;/*,j,cur;
	maxlen = (maxlen + 3) >>2;

	unsigned int ite;
	unsigned int curlen;

	char cap[6];
	cap[2] ='\0';
	uint32_t key;
	uint32_t nbsub;
	hashes.setSize(1);

	for(i=0;i< entries.size();i++){
		key =0;
		curlen = strlen(entries[i]);
		if (curlen > 0){
			cap[1] = entries[i][curlen-1];
			cap[3] = entries[i][0];
			if (curlen > 1){
				cap[0] = entries[i][curlen-2];
				cap[4] = entries[i][1];
				cap[5] = (curlen > 2) ? entries[i][2] : '\0';
			}else{
				cap[4] = '\0';
				cap[5] = '\0';
			}
		}else{
			cap[3] = '\0';
			cap[4] = '\0';
			cap[5] = '\0';
		}

		for(cur=0; entries[i][cur] != '\0';cur++){
			key =  *(uint32_t*)(curlen - cur > 2 ? entries[i] + cur : cap + 3 + cur - curlen);
			ite =hashes[0].find(key);
			if (ite == 0xFFFFFFFF){
				nbsub = hashes.getSize();
				hashes[0][key] = nbsub | 0x80000000;
				key = *(uint32_t*)(curlen - cur > 6 ? entries[i] + cur+4 : curlen > cur+4 ? cap + 7 + cur - curlen: cap+2);
				hashes.push_back();
				hashes[nbsub][key] = i;
			}else{
				nbsub = hashes[0].deref(ite);
				key = *(uint32_t*)(curlen - cur > 6 ? entries[i] + cur+4 : curlen > cur+4 ? cap + 7 + cur - curlen: cap+2);
				ite = hashes[nbsub].find(key);
				if (ite ==  0xFFFFFFFF) hashes[nbsub][key] = i;
				// else TODO!
			}
		}
	}*/

	//Compare_key ck;
	for(i=0;i<hashes.size();i++) hashes[i].keysort();
}



void Dictionary::findSubs(Vector<const char*> &fout, const char* query, uint32_t max_requested){
	Vector< uint32_t> foutindex;
	this->findIndexes(foutindex,query,max_requested);
	for(uint32_t i=0;i<foutindex.getSize();i++) fout.push_back(entries[foutindex[i]]);
}

void Dictionary::findIndexes(Vector<uint32_t> &fout, const char* query, uint32_t max_requested){
	uint32_t i = strlen(query);
	char* mquery;
	const char* c_start;
	const char* c_end;
	int c;
	if (flags & 1){
		mquery = new char[i+1];
		memcpy(mquery, query, i+1);
		for(c=0;mquery[c] != '\0';c++){
			if (((mquery[c] & 0xC0) == 0x40)&&(mquery[c] <= 'Z')) mquery[c] |= 0x20;
		}
		c_start = mquery;
		c_end = mquery;
	}else {
		c_start = query;
		c_end = query;

	}

	Vector< Tuple<uint32_t, 4u> > found_indexes;


	while(true){
		while(true){
			if ((*c_end == '\0')||(*c_end == ' ')||(*c_end == '\t')) break;
			c_end++;
		}
		if (c_end != c_start){
			found_indexes.push_back(this->findPrefix(c_start, (int)(c_end - c_start)));
		}
		if (*c_end == '\0') break;
		c_end++;
		c_start = c_end;
	}

	/*
	if (i < 3) hash_select =0;
	else{
		hash_select = hashes[0].find(*(uint32_t*)query);
		if (hash_select == -1) return;
		hash_select = hashes[0].deref(hash_select);
	}
	int j;
	if ((i & 3) == 0){
		for(j=0;j< hashes[hash_select].getSize();j++){
			if (hashes[hash_select].deref(j) & 0x80000000){
				// has subhash, todo
			}else{
				fout.push_back(entries[hashes[hash_select].deref(j)]);
				if (fout.size() == max_requested) return;
			}
		}
	}*/

	uint32_t j;
	int k,l;
	Vector< KeyElem<uint32_t, uint32_t> > scope;

	if (found_indexes.getSize() == 1){
		//printf("found indexes: %X %X %X %X\n", found_indexes[0][0], found_indexes[0][1], found_indexes[0][2],found_indexes[0][3]);
		if (found_indexes[0][0] != 0xFFFFFFFF)
		for(i=found_indexes[0][1];i<hashes[found_indexes[0][0]].getSize();i++){
			if ((hashes[found_indexes[0][0]].deref_key(i) ^ found_indexes[0][2]) & found_indexes[0][3]) break; // mismatch
			k = hashes[found_indexes[0][0]].deref(i);
			if (hashes[found_indexes[0][0]].deref_key(i) & 0xFF){ // sub is hash!
				scope.push_back( KeyElem<uint32_t, uint32_t>(k,0) );
				while( (l = scope.size()) > 0){
					l--;
					if (hashes[scope[l].k].deref_key(scope[l].d) & 0xFF){
						scope.push_back( KeyElem<uint32_t, uint32_t>(hashes[scope[l].k].deref(scope[l].d),0) );
					}else{
						k = hashes[scope[l].k].deref(scope[l].d);
						for(j=0 ;j <word_links[k].getSize();j++) {
							fout.push_back(word_links[k][j]);
							if (fout.getSize() == max_requested) {
								fout.sort_unique();
								if (fout.getSize() == max_requested){
									if (flags & 1) delete[](mquery);
									return;
								}
							}
						}
						while(true){
							scope[l].d++;
							if ( scope[l].d < hashes[scope[l].k].getSize()) break;
							else{
								scope.pop_back();
								if (l ==0) break;
								l--;
							}
						}
					}
				}
			}else{

				for(j=0 ;j <word_links[k].getSize();j++) {

					fout.push_back(word_links[k][j]);
					if (fout.getSize() == max_requested) {
						fout.sort_unique();
						if (fout.getSize() == max_requested){
							if (flags & 1) delete[](mquery);
							return;
						}
					}
				}
			}
		}
	}else{
		/*
		int *subind = new int[word_indexes.getSize()];
		uint32_t curmin;
		uint32_t curmin_index;

		for(i=0;i<word_indexes.getSize();i++) subind [i] =0;*/
	}

	if (fout.getSize() != 0) {
		//foutindex.show();
		fout.sort_unique();
		//printf("after sort uniq\n");
		//foutindex.show();
		if (flags & 1) delete[](mquery);
		return;
	}

	if (flags & 1) delete[](mquery);
	// found nothing, be clever ?
}
uint32_t Dictionary::findEntry(const char* query){
    // if query is a word, ez-pez;
    const char* c_start = query;
	const char* c_end = query;

    Vector<uint32_t> wordind;
    while(true){
		while(true){
			if ((*c_end == '\0')||(*c_end == ' ')||(*c_end == '\t')) break;
			c_end++;
		}
		if (c_end != c_start){
			wordind.push_back(this->findWord(c_start,(int)(c_end - c_start)));
			if (wordind.last() == 0xFFFFFFFF) return 0xFFFFFFFF;
		}
		if (*c_end == '\0') break;
		c_end++;
		c_start = c_end;
	}
    uint32_t i,j;
    for(i=0,j=1;j<wordind.getSize();j++){
        if (word_links[wordind[j]].getSize() < word_links[wordind[i]].getSize()) i=j;
    }
    j = wordind[i];
    for(i=0;i< word_links[j].getSize();i++) if (strcmp(query, entries[word_links[j][i]]) == 0) break;
    return (i< word_links[j].getSize()) ? word_links[j][i] : 0xFFFFFFFF;
}
Vector<uint32_t> Dictionary::findEntries(const vector<string> &queries){Vector<uint32_t> fout;
    fout.setSize(queries.size());
    uint32_t i,j;
    for(i=0,j=0;i< fout.getSize();i++)  if ((fout[i] = this->findEntry(queries[i+j].c_str())) == 0xFFFFFFFF) {printf("Warning, did not find %s\n", queries[i+j].c_str()); i--; j++; fout.pop_back();}
return fout;}

ERRCODE Dictionary::save(FILE*f) const{
    uint32_t i = entries.getSize();
    fwrite(&flags,sizeof(uint32_t),1,f);
    fwrite(&i,sizeof(uint32_t),1,f);
    for(i=0;i<entries.getSize();i++){
        fwrite(entries[i],sizeof(char),strlen(entries[i])+1,f);
    }
    word_links.save(f);
    return hashes.save(f);
}
ERRCODE Dictionary::load(FILE*f, unsigned int s){
    char buffer[256];
    char* cur;
    uint32_t i;
    fread(&flags,sizeof(uint32_t),1,f);
    fread(&i,sizeof(uint32_t),1,f);
    entries.setSize(i);
    for(i=0;i<entries.getSize();i++){
        cur = buffer;
        do{
            fread(cur,sizeof(char),1,f);
        }while (*(cur++) != '\0');
        entries[i] = new char[(int)(cur - buffer)];
        strcpy(entries[i], buffer);
    }
    word_links.load(f);
    hashes.load(f);
return 0;}

void Dictionary::show(FILE*f, int level) const{
	fprintf(f,"dictionary with %i entries\n", entries.size());
	for(uint32_t i=0;i<  entries.size();i++){
		fprintf(f,"%i: %s\n", i, entries[i]);
	}
	fprintf(f, " nb hashes %i \n",  hashes.getSize()  );

	/*uint32_t i;
	uint32_t j;
	uint32_t k;

	for(i=0;i< hashes[0].getSize();i++){
		j = *(uint32_t*)&hashes[0].deref_key(i);
		k = j;
		ExOp::bytereverse(k);
		fprintf(f,"[%.4s] [%X:%X] : %i entries\n", (char*)&hashes[0].deref_key(i),j,k, hashes[hashes[0].deref(i)].getSize());
	}*/

	//for(i=0;i< word_links.getSize();i++) fprintf(stdout, "%s: %i entries linked\n",   word_links[i].k, word_links[i].d.getSize() );
}



void Annotation::genDicos(Tuple<uint32_t> &fout){
    uint32_t i;
    for(i=0;i<fout.getSize();i++){
        fout[i] = dicos.getSize();
        dicos.push_back();
    }
}

void Annotation::useAnnotation()const{
    lfhstatic_ressources[(uint32_t) LFHSTATIC_RESSOURCE_ANNOTATIONS ] = (void*) this;
}

ERRCODE Annotation::read_Rtable(const char* path, TMatrix<double, 0u,0u>& target){
    FILE* f = fopen(path,"r+");
    if (f == NULL) return ERRCODE_NOFILE;
    char buffer[256];
    char bufferB[256];
    Vector<string> colnames;
    Vector<string> rownames;
    uint32_t i;
    do{
        if (2 != fscanf(f,"%[^\t\n ]%[\t\n\r ]", buffer, bufferB)) return 1;
        for(i=0; (bufferB[i] != '\0')&&(bufferB[i] != '\n');i++) if (i == 255) break;
        colnames.push_back(string(buffer));
    }while(bufferB[i] == '\0');
    Vector<double*> tempbuf;
    float wrbuf;//colnames.show();
    while(!feof(f)){
        if (2 != fscanf(f,"%[^\t\n ]%[\t\n\r ]", buffer, bufferB)) return 1;
        rownames.push_back(string(buffer));
        tempbuf.push_back(new double[colnames.getSize()]);
        for(uint32_t j =0; j < colnames.getSize();j++){
            if (2 != fscanf(f,"%f%[\t\n\r ]", &wrbuf, bufferB)) return 1;
            tempbuf.last()[j] = wrbuf;
        }
    } //rownames.show();
    target.setSizes(colnames.getSize(),rownames.getSize());
    Tuple<uint32_t, 2u> coor;
    for(coor[1] = 0; coor[1]< rownames.getSize();coor[1]++){
        for(coor[0] = 0; coor[0] < colnames.getSize();coor[0]++){
            target(coor) = tempbuf[coor[1]][coor[0]];
        }
        delete[](tempbuf[coor[1]]);
    }

    Tuple<uint32_t> indexes_input; indexes_input.setSize(2);
    indexes_input[0] = dicos.getSize(); indexes_input[1] = indexes_input[0]+1;
    dicos.push_back();
    for(i = 0 ; i< colnames.getSize();i++) dicos.last().addEntry(colnames[i].c_str());
    dicos.push_back();
    for(i = 0 ; i< rownames.getSize();i++) dicos.last().addEntry(rownames[i].c_str());
    axes.addEntry((void*)&target).toMemmove(indexes_input);
   //printf("at %X found ", (uint32_t)(void*)&target); axes[(void*)&target].show();
    return 0;
}
ERRCODE Annotation::write_Rtable(const char* path, const TMatrix<double, 0u,0u>& target) const{
    FILE* f = fopen(path,"w+");
    if (f == NULL) return ERRCODE_NOFILE;

    uint32_t cite = axes.find((void*)&target);
    if (cite == 0xFFFFFFFF) {printf("TMatrix has no annotation, nothing is written\n"); return(1);}

    for(uint32_t i =0; i < dicos[axes.deref(cite)[0]].entries.getSize(); i++) fprintf(f,"%s%c", dicos[axes.deref(cite)[0]].entries[i], i+1 == dicos[axes.deref(cite)[0]].entries.getSize() ? '\n' : '\t');
    Tuple<uint32_t, 2u> coor;
    for(coor[1] =0; coor[1] < dicos[axes.deref(cite)[1]].entries.getSize(); coor[1]++) {
        fprintf(f,"%s", dicos[axes.deref(cite)[1]].entries[coor[1]]);
        for(coor[0]=0;coor[0]<target.sizes[0];coor[0]++) fprintf(f,"\t%f", target(coor));
        fprintf(f,"\n");
    }
    fclose(f);
    return 0;
}


ERRCODE Annotation::read_Rtable(const char* path, DataGrid<double, 2u>& target){
    FILE* f = fopen(path,"r+");
    if (f == NULL) return ERRCODE_NOFILE;
    char buffer[256];
    char bufferB[256];
    Vector<string> colnames;
    Vector<string> rownames;
    uint32_t i;
    do{
        if (2 != fscanf(f,"%[^\t\n ]%[\t\n\r ]", buffer, bufferB)) return 1;
        for(i=0; (bufferB[i] != '\0')&&(bufferB[i] != '\n');i++) if (i == 255) break;
        colnames.push_back(string(buffer));
    }while(bufferB[i] == '\0');
    Vector<double*> tempbuf;
    float wrbuf;//colnames.show();
    while(!feof(f)){
        if (2 != fscanf(f,"%[^\t\n ]%[\t\n\r ]", buffer, bufferB)) return 1;
        rownames.push_back(string(buffer));
        tempbuf.push_back(new double[colnames.getSize()]);
        for(uint32_t j =0; j < colnames.getSize();j++){
            if (2 != fscanf(f,"%f%[\t\n\r ]", &wrbuf, bufferB)) return 1;
            tempbuf.last()[j] = wrbuf;
        }
    } //rownames.show();
    Tuple<uint32_t, 2u> coor; coor[0] = colnames.getSize(); coor[1] = rownames.getSize();
    target.setSizes(coor);
    for(coor[1] = 0; coor[1]< rownames.getSize();coor[1]++){
        for(coor[0] = 0; coor[0] < colnames.getSize();coor[0]++){
            target(coor) = tempbuf[coor[1]][coor[0]];
        }
        delete[](tempbuf[coor[1]]);
    }

    Tuple<uint32_t> indexes_input; indexes_input.setSize(2);
    indexes_input[0] = dicos.getSize(); indexes_input[1] = indexes_input[0]+1;
    dicos.push_back();
    for(i = 0 ; i< colnames.getSize();i++) dicos.last().addEntry(colnames[i].c_str());
    dicos.push_back();
    for(i = 0 ; i< rownames.getSize();i++) dicos.last().addEntry(rownames[i].c_str());
    axes.addEntry((void*)&target).toMemmove(indexes_input);
   //printf("at %X found ", (uint32_t)(void*)&target); axes[(void*)&target].show();
    return 0;
}
ERRCODE Annotation::write_Rtable(const char* path, const DataGrid<double, 2u>& target) const{
    FILE* f = fopen(path,"w+");
    if (f == NULL) return ERRCODE_NOFILE;

    uint32_t cite = axes.find((void*)&target);
    if (cite == 0xFFFFFFFF) {printf("s has no annotation, nothing is written\n"); return(1);}

    for(uint32_t i =0; i < dicos[axes.deref(cite)[0]].entries.getSize(); i++) fprintf(f,"%s%c", dicos[axes.deref(cite)[0]].entries[i], i+1 == dicos[axes.deref(cite)[0]].entries.getSize() ? '\n' : '\t');
    Tuple<uint32_t, 2u> coor;
    for(coor[1] =0; coor[1] < dicos[axes.deref(cite)[1]].entries.getSize(); coor[1]++) {
        fprintf(f,"%s", dicos[axes.deref(cite)[1]].entries[coor[1]]);
        for(coor[0]=0;coor[0]<target.dims[0];coor[0]++) fprintf(f,"\t%f", target(coor));
        fprintf(f,"\n");
    }
    fclose(f);
return 0;}

ERRCODE Annotation::save(FILE*f) const{ ERRCODE fout;
    fout = dicos.save(f);
return fout;}
ERRCODE Annotation::load(FILE*f, unsigned int s) {ERRCODE fout;
    fout = dicos.load(f);
return fout;}

#ifdef Rcpp_hpp
Rcpp::IntegerVector Metable::genFactor(uint32_t col) const{
}

Rcpp::IntegerVector Metable::genFactor(string col) const{
}
#endif


/** \brief Set cached density weights from parameters
 *
 * \return
 *
 */
void MarkovAttractor::populateCoefs(){
    const uint32_t nbd = transit_and_prior.sizes[0];
    coefs.setSize(nbd * nbd*2);
    double* currow;
    double rate;
    uint32_t k;


    for(uint32_t j=0;j<nbd;j++){
        currow[j * (nbd+1)*2|1] =0;
        for(k=0; k < nbd;k++){
            currow[j * (nbd+1)*2|1] += transit_and_prior(j,k);
        }
        currow[j * (nbd+1)*2|1] = transit_and_prior(j,j) / currow[j * (nbd+1)*2|1];
    }

    for(uint32_t j=0;j<nbd;j++){
        currow =  coefs() + (j *2 * nbd);
        for(uint32_t i=0;i<nbd;i++){
            if (i == j){
                currow[i << 1] = log(transit_and_prior(i,i));
                for(k=0; k < i;k++) currow[i << 1] -= transit_and_prior(i,k);
                for(k++; k < nbd;k++) currow[i << 1] -= transit_and_prior(i,k);
            }else{
                currow[i << 1] = currow[j * (nbd+1)*2|1] * transit_and_prior(j,i);//
                currow[(i << 1)|1] = transit_and_prior(j,i);//

            }
        }
    }
    coefs[1] =0;
    for(k=0; k < nbd;k++) coefs[1] += exp( currow[k * 2 * (nbd+1)]);
    coefs[1] = log(1.0 - coefs[1]); // probability of being transiant

}

Tuple<double> MarkovAttractor::computeLL(const SparseMatrix<double>& hidden)const{ Tuple<double> fout;
    fout.setSize(hidden.getNBcols());
    const double* tmpptr[2];
    const uint32_t nbd = transit_and_prior.sizes[0];
    for(uint32_t i=0;i< hidden.getNBcols();i++){
        switch(hidden.data[i].getSize()){
            case 1: fout[i] = coefs[hidden.data[i].deref_key(0) * ((nbd+1)*2) ]; break;
            case 2: tmpptr[0] = coefs() + (nbd + (hidden.data[i].deref_key(0) * nbd + hidden.data[i].deref_key(1)) * 2);
                tmpptr[1] = coefs() + (nbd + (hidden.data[i].deref_key(1) * nbd + hidden.data[i].deref_key(0)) * 2);
                fout[i] = coefs[1] + log(tmpptr[0][0] * exp(tmpptr[0][1] * hidden.data[i].deref(0)) + tmpptr[1][0] * exp(tmpptr[1][1] * hidden.data[i].deref(1)));
                break;
            default:
                ExCo<double>::toNegInfinity(fout[i]); // Too sad
        }
    }
}


MatrixDistribution::MatrixDistribution():stat_scope(NULL){}
void MatrixDistribution::setSizes(Tuple<unsigned int, 2u> coor, bool init_with_default_prior){
    rowprc.setSize(coor[0]);
    colprc.setSize(coor[1]);
    means.setSizes(coor[0],coor[1]);
    if (init_with_default_prior){
        rowprc.toOne();
        colprc.toOne();
        means.toZero();
    }
}
/*
void MatrixDistribution::setCovars(Trianglix<double>& rowcov, Trianglix<double>& colcov){
    Tuple<unsigned int, 2u> coor; coor[0] = rowcov.getSize(); coor[1] = colcov.getSize(); this->setSizes();
    rowprc.setSize(coor[0]);
    colprc.setSize(coor[1]);
}*/
double MatrixDistribution::EMHiddenMatrix(const SparseMatrix<double> &data, DirichletDistribution<0u>& rowprior, DirichletDistribution<0u>& colprior, DataGrid<double, 2u> &rowhidden, DataGrid<double, 2u> &colhidden){
    // NOT WORKING!
    Vector<uint32_t> row_size = data.compute_rowsizes();
    DataGrid<double, 2u> proj_X;
    Trianglix<double> row_w_cumul; row_w_cumul.setSize(rowprc.getSize()); // sum_i ZZ^t / Z^tUZ
 //   Trianglix<double> col_w_cumul; col_w_cumul.setSize(colprc.getSize()); // sum_j HH^t / H^tVH
    Tuple<double> meanproject;
    Tuple<double> rowsums[2];
    Tuple<double> colsums[2];

    Tuple<double> rowdenum[2];
    Tuple<double> coldenum[2];
    int cached_index;
    for(cached_index =0; cached_index <2;cached_index++){
        rowsums[cached_index].setSize(row_size.getSize());
        colsums[cached_index].setSize(data.getNBcols());
        rowdenum[cached_index].setSize(row_size.getSize());
        coldenum[cached_index].setSize(data.getNBcols());
    }cached_index =0;

    Tuple<uint32_t, 2u> coor, coor2, coorC, coorR;
    coor[0] = rowprc.getSize();
    coor[1] = colprc.getSize();
    meanproject.setSize(rowprc.getSize());
    proj_X.setSizes(coor);
    proj_X.toZero();


    ProgressBarPrint pbar(20);

    double tmp,tmpsum;

    Trianglix<double> dll_dU;
    Trianglix<double> dll_dV;
    TMatrix<double> dll_dM;


    DataGrid<double> dll_dZ_stat;
    DataGrid<double> dll_dH_stat;

    Tuple<double> dll_dA;
    Tuple<double> dll_dB;

    DomainMappedTuple<double,0u,LFHDOMAINMAP_POSITIVE_REAL> guess_A(rowprior.param);
    DomainMappedTuple<double,0u,LFHDOMAINMAP_POSITIVE_REAL> guess_B(colprior.param);


    dll_dU.setSize(rowhidden.dims[0]);
    dll_dV.setSize(colhidden.dims[0]);
    dll_dM.setSizes(rowhidden.dims[0], colhidden.dims[0]);
    coor[0] = colhidden.dims[0]; coor[1] = row_size.getSize(); dll_dZ_stat.setSizes(coor);
    coor[0] = rowhidden.dims[0]; coor[1] = data.getNBcols(); dll_dH_stat.setSizes(coor);
    dll_dA.setSize(rowhidden.dims[0]);
    dll_dB.setSize(colhidden.dims[0]);

    double ll;
    int i,k;

    double ll_best;

    Tuple<double,0u,TUPLE_FLAG_REMOTE_MEMORY> dainput;

    uint32_t rowite;
    int step;

    FuzzyCurvatureSearchScope<7> alpha;

    for(step=0;step<10;step++){
        // step 1: compute likelihood (only, but cache some calculation if pertinent for derivatives)
        ll = rowprior.getLL(rowhidden) + colprior.getLL(colhidden);

        printf("tmp_ll: %e\n", ll);
        dll_dZ_stat.toZero(); dll_dH_stat.toZero();
        dainput.setSize(rowhidden.dims[0]);
        for(coor[1]=0;coor[1]<rowhidden.dims[1];coor[1]++){
            coor[0]=0; dainput = rowhidden(coor); // set address to tuple
            if (coor[1] < 5) dainput.show();
           // for(coor[0]=0;coor[0]<rowhidden.dims[0];coor[0]++) dainput[coor[0] ] = rowhidden(coor);
            rowdenum[cached_index][coor[1]] = rowprc.Xformed_inner_product(dainput);
           // ll -= log(rowdenum[cached_index][coor[1]]) * data.getRowSize(coor[1]);
            //row_w_cumul += Trianglix<double>(dainput) / rowdenum[coor[1]];
        }
        dainput.setSize(colhidden.dims[0]);
        pbar.start("Compute LL");
        for(coor[1]=0;coor[1]<data.getNBcols();coor[1]++){ pbar.update(coor[1],data.getNBcols());// column
            coorC[1] = coor[1];
            coor[0]=0; dainput = colhidden(coor); // set address to tuple
            //for(coor[0]=0;coor[0]<colhidden.dims[0];coor[0]++) dainput[coor[0] ] = rowhidden(coor);
            coldenum[cached_index][coor[1]] = colprc.Xformed_inner_product(dainput);
          //  ll -= log(coldenum[cached_index][coor[1]]) * data.getColSize(coor[1]);
            if (!ExOp::isValid(log(coldenum[cached_index][coor[1]]))){

                printf("%e = %e\n", colprc.Xformed_inner_product(dainput), log(coldenum[cached_index][coor[1]]));
                colprc.show();
                dainput.show();
            }
            //col_w_cumul += Trianglix<double>(dainput) / coldenum;
            if ((coor[1] & 127) == 127) printf("tmp_ll: %e\n", ll);

            for(coorR[0]=0;coorR[0]<rowprc.getSize();coorR[0]++){
                coorC[0]=0;
                meanproject[ coorR[0] ] = means.data[coorR[0] ] * colhidden(coorC);
                for(coorC[0]++;coorC[0]<colprc.getSize();coorC[0]++) meanproject[ coorR[0]  ] +=  means.data[coorR[0]+ coorC[0] * rowprc.getSize()]  *colhidden(coorC);
            }


            for(rowite=0;rowite<data.data[coor[1]].getSize();rowite++){ // row
                coor[0] = data.data[coor[1]].deref_key(rowite);
                coorR[1] = coor[0];
                tmp = data.data[coor[1]].deref(rowite);
                for(coorR[0]=0;coorR[0]<rowprc.getSize();coorR[0]++){
                    tmp -= meanproject[coorR[0]] * rowhidden(coorR);
                }

                for(coorR[0]=0,coorC[0]=0;coorR[0] < colhidden.dims[0];coorR[0]++,coorC[0]++) dll_dZ_stat(coorR) += colhidden(coorC) * tmp / coldenum[cached_index][coorR[1]];
                for(coorR[0]=0,coorC[0]=0;coorC[0] < rowhidden.dims[0];coorR[0]++,coorC[0]++) dll_dH_stat(coorC) += rowhidden(coorR) * tmp / rowdenum[cached_index][coorC[1]];

                tmp = tmp * tmp;
                ll -= tmp / (rowdenum[cached_index][coor[0]] * coldenum[cached_index][coor[1]]);
                rowsums[cached_index][coor[0]] += tmp / rowdenum[cached_index][coor[0]];
                colsums[cached_index][coor[1]] += tmp / coldenum[cached_index][coor[1]];



                /*
                for(coor2[1]=0;coor2[1]<proj_X.dims[1];coor2[1]++){
                    coorC[0] = coor2[1];
                    tmp = data(coor) * colhidden(coorC);
                    for(coor2[0]=0;coor2[0]<proj_X.dims[0];coor2[0]++){
                        coorR[0] = coor2[0];
                        proj_X(coor2) += data(coor) * rowhidden(coorR);
                    }
                }*/
            }



        } pbar.finish();
        printf("current!\n");
        printf("LL[%i] = %e\n", step, ll);
        printf("P dLL_da = "); ExOp::show(dll_dA);
        tmp = this->evalLL_routine(data,rowprior,colprior,rowhidden,colhidden,means,rowprc,colprc);
        Tuple<double> testest; testest.setSize(dll_dA.getSize());
        for(coor[0]=0;coor[0]< dll_dA.getSize();coor[0]){
            rowprior.param[coor[0]] += 0.001f;
            rowprior.updateKonstant();
            testest[coor[0]] = this->evalLL_routine(data,rowprior,colprior,rowhidden,colhidden,means,rowprc,colprc);
            rowprior.param[coor[0]] -= 0.001f;
            testest[coor[0]] -= tmp;
        }
        testest *= 1000.0f;
        printf("R dLL_da = "); ExOp::show(testest);



        return 0.0f;
        if (alpha.needsUndo(ll, step > 0)){
            // undo all, BUT Z and H
            means -= dll_dM * (0.75f * alpha[0]);
            rowprc.toAddMult(dll_dU, -0.75f * alpha[1]);
            colprc.toAddMult(dll_dV, -0.75f * alpha[2]);


            guess_A.applyUpdate(dll_dA, -0.75f * alpha[5]);
            guess_B.applyUpdate(dll_dB, -0.75f * alpha[6]);
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


            rowprior.wrLLDerivative(rowhidden, dll_dA);
            colprior.wrLLDerivative(colhidden, dll_dB);

            dll_dU.toZero();
            for(coor[1]=0;coor[1]<rowhidden.dims[1];coor[1]++){
                for(i=0,k=0;i<dll_dU.getSize();i++){
                    coor[0] = i;
                    tmp = rowhidden(coor) * ( (rowsums[cached_index][i] / rowdenum[cached_index][coor[0]]) + data.getColSize(coor[1])) /rowdenum[cached_index][coor[0]];
                    for(coor[0]=0;coor[0]<=i;coor[0]++) dll_dU.data[k++] += rowhidden(coor) * tmp;
                }
            }
            dll_dV.toZero();
            for(coor[1]=0;coor[1]<colhidden.dims[1];coor[1]++){
                for(i=0,k=0;i<dll_dV.getSize();i++){
                    coor[0] = i;
                    tmp = colhidden(coor) * ( (colsums[cached_index][i] / coldenum[cached_index][coor[0]]) + row_size[coor[1]] ) /coldenum[cached_index][coor[0]];
                    for(coor[0]=0;coor[0]<=i;coor[0]++) dll_dV.data[k++] += colhidden(coor) * tmp;
                }
            }

            cached_index = 1 - cached_index;
            means += dll_dM * alpha[0];
            rowprc.toAddMult(dll_dU, alpha[1]);
            colprc.toAddMult(dll_dV, alpha[2]);

            guess_A.applyUpdate(dll_dA, alpha[5]);
            guess_B.applyUpdate(dll_dB, alpha[6]);
        }
        rowprior.updateKonstant();colprior.updateKonstant();
    }while(false);
return 0.0;}
double MatrixDistribution::EMHiddenMatrix(const SparseMatrix<double> &data, DataGrid<double, 2u> &rowhidden, DataGrid<double, 2u> &colhidden){
    DataGrid<double, 2u> proj_X;
    Vector<uint32_t> row_size = data.compute_rowsizes();
    Trianglix<double> row_w_cumul; row_w_cumul.setSize(rowprc.getSize()); // sum_i ZZ^t / Z^tUZ
 //   Trianglix<double> col_w_cumul; col_w_cumul.setSize(colprc.getSize()); // sum_j HH^t / H^tVH
    Tuple<double> meanproject;
    Tuple<double> rowsums[2];
    Tuple<double> colsums[2];
    TMatrix<double> Msums[2];
    Tuple<double> Msumstmp;
    Tuple<double> rowdenum[2];
    Tuple<double> coldenum[2];

    DataGrid<double, 2u> rowh[2];
    DataGrid<double, 2u> colh[2];

    int cached_index;
    for(cached_index =0; cached_index <2;cached_index++){

        rowsums[cached_index].setSize(row_size.getSize());
        colsums[cached_index].setSize(data.getNBcols());
        rowdenum[cached_index].setSize(row_size.getSize());
        coldenum[cached_index].setSize(data.getNBcols());
        Msums[cached_index].setSizes(rowhidden.dims[0],colhidden.dims[0]);
    }cached_index =0;
    rowh[1] = rowhidden;
    rowh[0].setSizes(rowh[1].dims);
    colh[1] = colhidden;
    colh[0].setSizes(colh[1].dims);

    Msumstmp.setSize(rowhidden.dims[0]);

    Tuple<uint32_t, 2u> coor, coor2, coorC, coorR;
    coor[0] = rowprc.getSize();
    coor[1] = colprc.getSize();
    meanproject.setSize(rowprc.getSize());
    proj_X.setSizes(coor);
    proj_X.toZero();


    ProgressBarPrint pbar(20);

    double tmp,tmpsum;

    Trianglix<double> dll_dU;
    Trianglix<double> dll_dV;


    DataGrid<double> dll_dZ_stat;
    //DataGrid<double> dll_dH_stat;




    dll_dU.setSize(rowhidden.dims[0]);
    dll_dV.setSize(colhidden.dims[0]);

    coor[0] = colhidden.dims[0]; coor[1] =row_size.getSize(); dll_dZ_stat.setSizes(coor);
   // coor[0] = rowhidden.dims[0]; coor[1] = data.getNBcols(); dll_dH_stat.setSizes(coor);

    double ll;

    int i,k;

    double ll_best;

    Tuple<double,0u,TUPLE_FLAG_REMOTE_MEMORY> dainput;

    uint32_t rowite;
    int step;

    FuzzyCurvatureSearchScope<3> alpha;

    for(step=0;step<10;step++){
        // step 1: compute likelihood (only, but cache some calculation if pertinent for derivatives)
        ll = 0.0;
        dll_dZ_stat.toZero(); // dll_dH_stat.toZero();
        dainput.setSize(rowhidden.dims[0]);
        rowsums[cached_index].toZero();
        colsums[cached_index].toZero();
        Msums[cached_index].toZero();
        for(coor[1]=0;coor[1]<rowhidden.dims[1];coor[1]++){
            if (row_size[coor[1]] == 0) continue; // nothing to see...
            coor[0]=0; dainput = rowh[cached_index^1](coor); // set address to tuple
           // for(coor[0]=0;coor[0]<rowhidden.dims[0];coor[0]++) dainput[coor[0] ] = rowhidden(coor);
            rowdenum[cached_index][coor[1]] = rowprc.Xformed_inner_product(dainput);
            ll -= log(rowdenum[cached_index][coor[1]]) * row_size[coor[1]];
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
            /*if (!ExOp::isValid(log(coldenum[cached_index][coor[1]]))){
                printf("%e = %e\n", colprc.Xformed_inner_product(dainput), log(coldenum[cached_index][coor[1]]));
                colprc.show();
                dainput.show();
                LFH_ALIVE; exit(1);
            }*/
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
            rowprc.toAddMult(dll_dU, -0.75f * alpha[1]);
            colprc.toAddMult(dll_dV, -0.75f * alpha[2]);
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
                tmpsum = ( (rowsums[cached_index][coor[1]] / rowdenum[cached_index][coor[1]]) - row_size[coor[1]]) / rowdenum[cached_index][coor[1]];
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

            /*printf("P dLL_dU = "); ExOp::show(dll_dU);
            printf("P dLL_dV = "); ExOp::show(dll_dV);
            printf("P dLL_dM = "); ExOp::show(Msums[cached_index]);
            tmp = this->evalLL_routine(data,rowhidden,colhidden,means,rowprc,colprc);
            Trianglix<double> testest; testest.setSize(dll_dU.getSize());
            printf("reeval LL = %e\n", tmp);*/
            /*for(coor[0]=0,i=0;coor[0]< dll_dU.getSize() ;coor[0]++){
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
            */
            /*TMatrix<double> testest2; testest2.setSizes(rowhidden.dims[0],colhidden.dims[0]);
            for(coor[0]=0,i=0;coor[0]< colhidden.dims[0] ;coor[0]++){
                for(coor[1]=0; coor[1]< rowhidden.dims[0]; coor[1]++,i++){
                    means.data[i] += 0.000001f;
                    testest2.data[i] = this->evalLL_routine(data,rowhidden,colhidden,means,rowprc,colprc) - tmp;
                    means.data[i] -= 0.000001f;
                }
            }
            testest2 *= 1000000.0f;
            printf("R dLL_dM = "); ExOp::show(testest2);
            return 0.0f;*/


            means +=  Msums[cached_index] * alpha[0];
            rowprc.toAddMult(dll_dU, alpha[1]);
            colprc.toAddMult(dll_dV, alpha[2]);
            cached_index = 1 - cached_index;
        }
        break;
    }
return 0.0;}

double MatrixDistribution::evalLL_routine(const SparseMatrix<double> &data, DirichletDistribution<0u>& rowprior, DirichletDistribution<0u>& colprior, DataGrid<double, 2u> &rowhidden, DataGrid<double, 2u> &colhidden, const TMatrix<double,0u,0u> &mat_M, const Trianglix<double> &mat_U, const Trianglix<double> &mat_V)const{
    Tuple<unsigned int, 2u> coor, coorC, coorR;
    double fout = 0.0f; // rowprior.getLL(rowhidden) + colprior.getLL(colhidden);
    uint32_t rowite;
    Tuple<double,0u,TUPLE_FLAG_REMOTE_MEMORY> dainput;
    dainput.setSize(rowhidden.dims[0]);
    Vector<uint32_t> rowsize = data.compute_rowsizes();
    Tuple<double> rowdenum; rowdenum.setSize(rowsize.getSize());
    Tuple<double> coldenum; coldenum.setSize(data.getNBcols());
    Tuple<double> meanproject;
    double tmp;
    meanproject.setSize(rowprc.getSize());
    for(coor[1]=0;coor[1]<rowsize.getSize();coor[1]++){
        coor[0]=0; dainput = rowhidden(coor);
        rowdenum[coor[1]] = mat_U.Xformed_inner_product(dainput);
        //fout -= log(rowdenum[coor[1]]) * data.getRowSize(coor[1]);
    }
    dainput.setSize(colhidden.dims[0]);
    for(coor[1]=0;coor[1]<data.getNBcols();coor[1]++){ // column
        coor[0]=0; dainput = colhidden(coor); // set address to tuple
        coldenum[coor[1]] = mat_V.Xformed_inner_product(dainput);
        //fout -= log(coldenum[coor[1]]) * data.getColSize(coor[1]);

        for(coorR[0]=0;coorR[0]<rowprc.getSize();coorR[0]++){
            coorC[0]=0;
            meanproject[ coorR[0] ] = means.data[coorR[0] ] * colhidden(coorC);
            for(coorC[0]++;coorC[0]<colprc.getSize();coorC[0]++) meanproject[ coorR[0]  ] += mat_M.data[coorR[0]+ coorC[0] * rowprc.getSize()]  *colhidden(coorC);
        }

        for(rowite=0;rowite<data.data[coor[1]].getSize();rowite++){ // row
            coorR[1] = coor[0] = data.data[coor[1]].deref_key(rowite);
            tmp = data.data[coor[1]].deref(rowite);
            for(coorR[0]=0;coorR[0]<rowprc.getSize();coorR[0]++){
                tmp -= meanproject[coorR[0]] * rowhidden(coorR);
            }
            tmp = 0.5 * tmp * tmp;
            fout -= tmp / (rowdenum[coor[0]] * coldenum[coor[1]]);
        }
    }
return fout;}
double MatrixDistribution::evalLL_routine(const SparseMatrix<double> &data, DataGrid<double, 2u> &rowhidden, DataGrid<double, 2u> &colhidden, const TMatrix<double,0u,0u> &mat_M, const Trianglix<double> &mat_U, const Trianglix<double> &mat_V)const{
    Tuple<unsigned int, 2u> coor, coorC, coorR;
    double fout = 0.0;
    uint32_t rowite;
    Vector<uint32_t> rowsize = data.compute_rowsizes();

    Tuple<double,0u,TUPLE_FLAG_REMOTE_MEMORY> dainput;
    dainput.setSize(rowhidden.dims[0]);
    Tuple<double> rowdenum; rowdenum.setSize(rowsize.getSize());
    Tuple<double> coldenum; coldenum.setSize(data.getNBcols());
    Tuple<double> meanproject;
    double tmp;
    meanproject.setSize(rowhidden.dims[0]);
    for(coor[1]=0;coor[1]<rowsize.getSize();coor[1]++){
        coor[0]=0; dainput = rowhidden(coor);
        rowdenum[coor[1]] = mat_U.Xformed_inner_product(dainput);
        fout -= log(rowdenum[coor[1]]) * rowsize[coor[1]];
    }
    dainput.setSize(colhidden.dims[0]);
    for(coor[1]=0;coor[1]<data.getNBcols();coor[1]++){ // column
        coorC[1] = coor[1];
        coor[0]=0; dainput = colhidden(coor); // set address to tuple
        coldenum[coor[1]] = mat_V.Xformed_inner_product(dainput);
        fout -= log(coldenum[coor[1]]) * data.getColSize(coor[1]);

        for(coorR[0]=0;coorR[0]<mat_U.getSize();coorR[0]++){
            coorC[0]=0;
            meanproject[ coorR[0] ] = means.data[coorR[0] ] * colhidden(coorC);
            for(coorC[0]++;coorC[0]<colprc.getSize();coorC[0]++) meanproject[ coorR[0]  ] += mat_M.data[coorR[0]+ coorC[0] * rowprc.getSize()]  *colhidden(coorC);
        }

        for(rowite=0;rowite<data.data[coor[1]].getSize();rowite++){ // row
            coor[0] = data.data[coor[1]].deref_key(rowite);
            tmp = data.data[coor[1]].deref(rowite);
            coorR[1] = coor[0];
            for(coorR[0]=0;coorR[0]<mat_U.getSize();coorR[0]++){
                tmp -= meanproject[coorR[0]] * rowhidden(coorR);
            }
            tmp = 0.5 * tmp * tmp;
            //if (rowite < 3) printf("%i: %e\n", rowite, tmp);
            fout -= tmp / (rowdenum[coor[0]] * coldenum[coor[1]]);
        }
    }
return fout;}

double MatrixDistribution::LL(const DataGrid<double, 2u>& value) const{
    double fout;
    unsigned int i,j,k;

    double* buffer = new double[value.dims[0]*value.dims[1]];

    Tuple<unsigned int, 2u> coor;
    double* minibuf = new double[value.dims[0] < value.dims[1] ? value.dims[1] : value.dims[0]];
    for(coor[1]=0;coor[1]<value.dims[1];coor[1]++){
        for(coor[0]=0;coor[0] < value.dims[0];coor[0]++) minibuf[coor[0]] = value(coor) - means.data[coor[1] * rowprc.getSize() +coor[0]];
        for(i=0,k=0;i<rowprc.getSize();i++,k++){
            buffer[coor[1] * value.dims[0] +i] =0;
            for(j=0;j<i;j++,k++){
                buffer[coor[1] * value.dims[0] +i] += rowprc.data[k] * minibuf[j];
                buffer[coor[1] * value.dims[0] +j] += rowprc.data[k] * minibuf[i];
            }
            buffer[coor[1] * value.dims[0] +i] += rowprc.data[k] * minibuf[i];
        }
    }
    fout =0.0f;
    for(coor[0]=0;coor[0] < value.dims[0];coor[0]++){
        for(i=0,k=0;i<colprc.getSize();i++,k++){
            minibuf[i] =0;
            for(j=0;j<i;j++,k++){
                minibuf[i] += colprc.data[k] * buffer[j * value.dims[0] + coor[0]];
                minibuf[j] += colprc.data[k] * buffer[i * value.dims[0] + coor[0]];
            }
            minibuf[i] += colprc.data[k] * buffer[i * value.dims[0] + coor[0]];
        }
        for(i=0;i<colprc.getSize();i++) fout -= minibuf[i] * minibuf[i];
    }
    delete[](buffer);
    delete[](minibuf);

return fout;}
void MatrixDistribution::EMinit(){
    if (stat_scope == NULL){
        stat_scope = new KeyElem< Tuple<Trianglix<double>, 2u>, double* >();
    }else delete[](stat_scope->d);
    stat_scope->d = new double[(rowprc.getSize() +1)* (colprc.getSize()+1)];
    stat_scope->k[0].setSize(rowprc.getSize());
    stat_scope->k[1].setSize(colprc.getSize());
    ExOp::toZero(stat_scope->k);
    for(uint32_t i=0;i<(rowprc.getSize() +1)* (colprc.getSize()+1);i++) stat_scope->d[i] = 0.0f;
}
void MatrixDistribution::EMregist(const DataGrid<double, 2u> &instance, const double prob){

}

double MatrixDistribution::EMHiddenMatrix(const DataGrid<double, 2u> &data, const DataGrid<double, 2u> &rowhidden, const DataGrid<double, 2u> &colhidden){
    DataGrid<double, 2u> proj_X;
    Trianglix<double> row_w_cumul; row_w_cumul.setSize(rowprc.getSize());
    Trianglix<double> col_w_cumul; col_w_cumul.setSize(colprc.getSize());
    Tuple<double> rowsums,colsums,meanproject; rowsums.setSize(rowprc.getSize()); colsums.setSize(colprc.getSize());
    Tuple<uint32_t, 2u> coor, coor2, coorC, coorR;
    coor[0] = rowprc.getSize();
    coor[1] = colprc.getSize();
    meanproject.setSize(rowprc.getSize());
    proj_X.setSizes(coor);
    proj_X.toZero();

    double tmp;

    Tuple<double> rowdenum; rowdenum.setSize(rowprc.getSize());
    double coldenum;

    double logdet_row;
    double logdet_col;
    double ll;

    Tuple<double> dainput;


    do{
        logdet_row = 0.0;
        logdet_col = 0.0;
        ll = 0.0;
        dainput.setSize(rowhidden.dims[0]);
        for(coor[1]=0;coor[1]<rowhidden.dims[1];coor[1]++){
            for(coor[0]=0;coor[0]<rowhidden.dims[0];coor[0]++){
                dainput[rowhidden.dims[0] ] = rowhidden(coor);
            }
            rowdenum[coor[1]] = rowprc.Xformed_inner_product(dainput);
            logdet_row -= log(rowdenum[coor[1]]);
            row_w_cumul += Trianglix<double>(dainput) / rowdenum[coor[1]];
        }

        dainput.setSize(colhidden.dims[0]);

        for(coor[1]=0;coor[1]<data.dims[1];coor[1]++){ // column
            coorC[1] = coor[1];
            for(coor[0]=0;coor[0]<colhidden.dims[0];coor[0]++){
                dainput[colhidden.dims[0] ] = colhidden(coor);
            }
            coldenum = colprc.Xformed_inner_product(dainput);
            logdet_col -= log(coldenum);
            col_w_cumul += Trianglix<double>(dainput) / coldenum;
            for(coorR[0]=0;coorR[0]<rowprc.getSize();coorR[0]++){
                coorC[0]=0;
                meanproject[ coorR[0] ] = means.data[coorR[0] ] *colhidden(coorC);
                for(;coorC[0]<colprc.getSize();coorC[0]++) meanproject[ coorR[0]  ] +=  means.data[coor2[0]+ coorC[0]* rowprc.getSize()]  *colhidden(coorC);
            }
            for(coor[0]=0;coor[0]<data.dims[0];coor[0]++){ // row
                coorR[1] = coor[0];
                tmp = data(coor);
                for(coorR[0]=0;coorR[0]<proj_X.dims[0];coorR[0]++){
                    tmp -= meanproject[coorR[0]] * rowhidden(coorR);
                }
                tmp = tmp * tmp;
                ll -= tmp / (rowdenum[coor[0]] * coldenum);
                rowsums[coor[0]] += tmp / rowdenum[coor[0]];
                colsums[coor[1]] += tmp / coldenum;

                for(coor2[1]=0;coor2[1]<proj_X.dims[1];coor2[1]++){
                    coorC[0] = coor2[1];
                    tmp = data(coor) * colhidden(coorC);
                    for(coor2[0]=0;coor2[0]<proj_X.dims[0];coor2[0]++){
                        coorR[0] = coor2[0];
                        proj_X(coor2) += data(coor) * rowhidden(coorR);
                    }
                }
            }
        }
        ll += logdet_row * data.dims[1] + logdet_col * data.dims[0];
        printf("Log likelihood: %e\n", ll);
    }while(false);
return 0.0;}
/*
void MatrixDistribution::EMregistOuter(const double *rowhidden , const double* &colhidden, double value, const double prob){
    uint32_t i,j;
    double tmp;
    for(j=0;j<colprc.getSize();j++){
        tmp = value * rowhidden[j] * prob;
        for(i=0;i<rowprc.getSize();i++) stat_scope->d[i + j * rowprc.getSize()] += tmp *colhidden[i];
    }
    for(i=0;i<rowprc.getSize();i++) stat_scope->d[i + j * rowprc.getSize()] += rowhidden[i] * value;
    j++;
    for(i=0;i<colprc.getSize();i++) stat_scope->d[i + j * rowprc.getSize()] += colhidden[i] * value;
    stat_scope->d[i + j * rowprc.getSize()] += value;
}*/
double MatrixDistribution::EMfinit(){

    return 0.0;
}

void EigenFixer::operator()(Tuple<double> & eigenvals) const{
    // makes eigenvalues strictly positive, models positive eigen vals as chi-square di
    double logeigen,tmp;
    uint32_t nbval=0;
    for(uint32_t i=0;i <eigenvals.getSize();i++){
        tmp = log(fabs(eigenvals[i]));
        if (ExOp::isValid(tmp)){nbval++;logeigen += tmp;}
    }
    if (nbval == 0) eigenvals.toOne();
    else{ logeigen  /= nbval;
        logeigen = exp(logeigen - 16.0); // machine epsilon * geomean eigen
        for(uint32_t i=0;i <eigenvals.getSize();i++){
            eigenvals[i] = fabs(eigenvals[i]);
            if (eigenvals[i] < logeigen) eigenvals[i] = logeigen;
        }
    }
}
/*  if ((i >0)&&(value > old[0])|| (!ExOp::isValid(value))){
            printf("F[%e] = %e Fail\n", modeseek, value);
            next = (modeseek * 7 + old[1]) / 8.0;
            modeseek = old[1] + 0.00000005f;
            double deriv = lngamma(1.0 + modeseek) + lngamma(1.0 + n - modeseek) + lngamma(1.0 + K - modeseek) + lngamma(1.0 + modeseek + N - K - n) - lncumNormal((modeseek - dev) / sigma);
            double d2eriv = d_lngamma_dx(1.0 + modeseek) - d_lngamma_dx(1.0 + n - modeseek) - d_lngamma_dx(1.0 + K - modeseek) + d_lngamma_dx(1.0 + modeseek + N - K - n) - d_lncumNormal_dx((modeseek - dev) / sigma) / sigma;
            modeseek = old[1] - 0.00000005f;
            deriv -= lngamma(1.0 + modeseek) + lngamma(1.0 + n - modeseek) + lngamma(1.0 + K - modeseek) + lngamma(1.0 + modeseek + N - K - n) - lncumNormal((modeseek - dev) / sigma);
            d2eriv -= d_lngamma_dx(1.0 + modeseek) - d_lngamma_dx(1.0 + n - modeseek) - d_lngamma_dx(1.0 + K - modeseek) + d_lngamma_dx(1.0 + modeseek + N - K - n) - d_lncumNormal_dx((modeseek - dev) / sigma) / sigma;
            modeseek = old[1];
            eval2[0] = d_lngamma_dx(1.0 + modeseek) - d_lngamma_dx(1.0 + n - modeseek) - d_lngamma_dx(1.0 + K - modeseek) + d_lngamma_dx(1.0 + modeseek + N - K - n) - d_lncumNormal_dx((modeseek - dev) / sigma) / sigma;
            eval2[1] = d2_lngamma_dx2(1.0 + modeseek) + d2_lngamma_dx2(1.0 + n - modeseek) + d2_lngamma_dx2(1.0 + K - modeseek) + d2_lngamma_dx2(1.0 + modeseek + N - K - n) - d2_lncumNormal_dx2((modeseek - dev) / sigma) / (sigma * sigma);
            //d_lngamma_dx(1.0 + modeseek) - d_lngamma_dx(1.0 + n - modeseek) - d_lngamma_dx(1.0 + K - modeseek) + d_lngamma_dx(1.0 + modeseek + N - K - n) - d_lncumNormal_dx((modeseek - dev) / sigma) / sigma;
            printf(" deriv %e vs %e\n", eval2[0] , deriv * 10000000.0);
            printf("d2eriv %e vs %e\n", eval2[1] , d2eriv * 10000000.0);
            modeseek = old[1] + 0.00005f;
            deriv = - lncumNormal((modeseek - dev) / sigma);
            d2eriv = - d_lncumNormal_dx((modeseek - dev) / sigma) / sigma;
            modeseek = old[1] - 0.00005f;
            deriv -=  - lncumNormal((modeseek - dev) / sigma);
            d2eriv -=  - d_lncumNormal_dx((modeseek - dev) / sigma) / sigma;
            modeseek = old[1];
            eval2[0] = - d_lncumNormal_dx((modeseek - dev) / sigma) / sigma;
            eval2[1] = - d2_lncumNormal_dx2((modeseek - dev) / sigma) / (sigma * sigma);
doublechk = eval2[0];
            printf(" deriv %e vs %e\n", eval2[0] , deriv * 10000.0);
            printf("d2eriv %e vs %e\n", eval2[1] , d2eriv * 10000.0);
            modeseek = old[1] + 0.00005f;
            deriv = lngamma(1.0 + modeseek) + lngamma(1.0 + modeseek + N - K - n);
            d2eriv = d_lngamma_dx(1.0 + modeseek) + d_lngamma_dx(1.0 + modeseek + N - K - n);
            modeseek = old[1] - 0.00005f;
            deriv -= lngamma(1.0 + modeseek) + lngamma(1.0 + modeseek + N - K - n);
            d2eriv -= d_lngamma_dx(1.0 + modeseek) + d_lngamma_dx(1.0 + modeseek + N - K - n);
            modeseek = old[1];
            eval2[0] = d_lngamma_dx(1.0 + modeseek) + d_lngamma_dx(1.0 + modeseek + N - K - n) ;
            eval2[1] = d2_lngamma_dx2(1.0 + modeseek) + d2_lngamma_dx2(1.0 + modeseek + N - K - n);
doublechk += eval2[0];
            printf(" deriv %e vs %e\n", eval2[0] , deriv * 10000.0);
            printf("d2eriv %e vs %e\n", eval2[1] , d2eriv * 10000.0);

            modeseek = old[1] + 0.00005f;
            deriv = lngamma(1.0 + n - modeseek) + lngamma(1.0 + K - modeseek);
            d2eriv = -d_lngamma_dx(1.0 + n - modeseek) - d_lngamma_dx(1.0 + K - modeseek);
            modeseek = old[1] - 0.00005f;
            deriv -= lngamma(1.0 + n - modeseek) + lngamma(1.0 + K - modeseek);
            d2eriv -= -d_lngamma_dx(1.0 + n - modeseek) - d_lngamma_dx(1.0 + K - modeseek);
            modeseek = old[1];
            eval2[0] = -d_lngamma_dx(1.0 + n - modeseek) - d_lngamma_dx(1.0 + K - modeseek);
            eval2[1] = d2_lngamma_dx2(1.0 + n - modeseek) + d2_lngamma_dx2(1.0 + K - modeseek);
doublechk += eval2[0];
            printf(" deriv %e vs %e\n", eval2[0] , deriv * 10000.0);
            printf("d2eriv %e vs %e\n", eval2[1] , d2eriv * 10000.0);
            modeseek = next;
            printf("double check : %e\n", doublechk);
        }else{*/
double HypergeomericNormalLogitPval(double dev, double sigma, uint32_t n, uint32_t K, uint32_t N, bool show){ // logit[ 0.5 + 0.5*\sum_k=0^K P(k) * erf( (dev -k) / (sigma *sqrt(2))) ]
        // to find the modd suppose the hypergeo is a normal ditribution, its mean and var is:
        // mu = nK / N;
        // var = mu * (N-K) * (N-n) / (N * (N+1))
        // hence this approximation becomes and integral: \int_k=-inf^K=inf exp(-(k - mu)^2/(2var^2)) * erf( (dev -k) / sigma)
        // this is equivalent to erf( (dev -k) / (sigma * sqrt(1 + 2*dev/sigma ) )

        // first find the sign of Logit pval...


        // first, find the maximal term...
    double modeseek = (((double)(n + 1)) *(K +1)) / (N +2);
    bool flip = ((modeseek < dev)^(sigma > 0.0));
    bool outside;

    if (sigma == 0.0) {
        if (dev < 0) dev = 0;
        else{
            if (dev > n) dev = n;
            if (dev > K) dev = K;
        }
        if (dev + N < n + K) dev = (n + K) - N;
        modeseek = -HypergeomericLogitPval(dev,n,K,N);
        if (!ExOp::isValid(modeseek)) printf("direct Hypergeo failed with args %e %i %i %i\n", dev,n,K,N);
        return(modeseek);
    }
    if (flip) sigma = -sigma;
    if (show) printf("guess mode flippy?: %e \n", modeseek); // 1 - p = sum log(k_i 1) - e_i
    double old[2]; double next; double doublechk;
    double eval2[2];

    class LocalFunctor{
        public:
        double n,K,N,dev,sigma;
        double operator()(const double& k) const{
            return d_lngamma_dx(1.0 + k) - d_lngamma_dx(1.0 + n - k) - d_lngamma_dx(1.0 + K - k) + d_lngamma_dx(1.0 + k + N - K - n) - d_lncumNormal_dx((k - dev) / sigma) / sigma;
        }
    };

    LocalFunctor derivfunc;
    derivfunc.n = n;
    derivfunc.K = K;
    derivfunc.N = N;
    derivfunc.dev = dev;
    derivfunc.sigma = sigma;

    old[0] = dev;
    if (old[0]< 0.0) old[0] = 0.0;
    if (old[0] + N < K + n) old[0] = K + n - N;
    if (old[0] > n) old[0] = n;
    if (old[0] > K) old[0] = K;
    modeseek  = monotonicSolver<LocalFunctor, double>(derivfunc, old[0], modeseek);

    /*
    for(unsigned int i=0;i<10;i++){
        double value = lngamma(1.0 + modeseek) + lngamma(1.0 + n - modeseek) + lngamma(1.0 + K - modeseek) + lngamma(1.0 + modeseek + N - K - n) - lncumNormal((modeseek - dev) / sigma);
        old[1] = modeseek;
        old[0] = value;
        double deriv = derivfunc(modeseek);
        double d2eriv = d2_lngamma_dx2(1.0 + modeseek) + d2_lngamma_dx2(1.0 + n - modeseek) + d2_lngamma_dx2(1.0 + K - modeseek) + d2_lngamma_dx2(1.0 + modeseek + N - K - n) - d2_lncumNormal_dx2((modeseek - dev) / sigma) / (sigma * sigma);
        modeseek -= deriv / d2eriv;
    }*/

    modeseek = floor(modeseek);

    Magscaled<double> summy;
    Magscaled<double> cur;

    // find which is the lead term

    cumNormal(summy, (modeseek - dev) / sigma);
    cumNormal(cur, (modeseek + 1.0 - dev) / sigma);

    double tmpsum[2];
    tmpsum[0] = log(summy.value) + summy.exponent + lngammachoose(modeseek, K) + lngammachoose(n - modeseek, N - K) -  lngammachoose(n, N);
    tmpsum[1] = log(cur.value) + cur.exponent + lngammachoose(modeseek+1.0, K) + lngammachoose(n - modeseek - 1.0, N - K) -  lngammachoose(n, N);

    //if (show){//
    if (log(cur.value) - log(summy.value) + cur.exponent - summy.exponent - log(modeseek+1) + log(((double)K)-modeseek) + log(((double)n) - modeseek) - log(((double)N-K-n)+ modeseek+1) > 0){
        modeseek = modeseek + 1.0;
        summy = cur;
    }
    //if ((diff > 0 )^(tmpsum[1] > tmpsum[0])){ // log(summy.value) - log(cur.value) + summy.exponent - cur.exponent
    //    printf("%e is not %e  (%e,%e,%e,%e)\n", diff, tmpsum[1] - tmpsum[0], log(modeseek+1), log(K-modeseek), log(n - modeseek), log(((double)N-K-n)+ modeseek+1));
    //    printf("1(%e vs %e)\n",  + log(modeseek+1) - log(((double)K)-modeseek)  - log(((double)n) - modeseek) + log(((double)N-K-n)+ modeseek+1),  log(summy.value) - log(cur.value) + summy.exponent - cur.exponent);
    //}
    //}

    if (show){
        for(uint32_t o =0; o<5;o++){
        cumNormal(cur, (modeseek + 2 - o - dev) / sigma);
        printf("F[%i][%e]  = %e\n", (uint32_t)(modeseek + 2 - o) ,(modeseek + 2 - o - dev) / sigma, log(cur.value) + cur.exponent + lngammachoose(modeseek + 2 - o, K) + lngammachoose(n - modeseek + 2 - o, N - K));
    }   }
   // printf("%i, summy %e    ", (int)modeseek, log_e_x_m1(-log(summy.value) -summy.exponent) ); summy.show();

    if (show) printf("init term: cnormal(%e) = %e * e^%e for mode %i\n", (modeseek - dev) / sigma, summy.value, summy.exponent, (int)modeseek);
    double hyperterm[2]; hyperterm[0] = 1.0; hyperterm[1] = 1.0;
    double dblk[2]; dblk[0] = modeseek; dblk[1] = modeseek;
    old[1] = summy.exponent;
    uint32_t i=0;
    do{
        if ((dblk[0] > 0.99)&&(dblk[0] + N - K - n > 0.99)){
            hyperterm[0] *= dblk[0];
            hyperterm[0] *= dblk[0] + N - K - n;
            dblk[0]-= 1.0;
            hyperterm[0] /= -dblk[0] + K;
            hyperterm[0] /= -dblk[0] + n;
            cumNormal(cur, (dblk[0] - dev) / sigma);
            cur *= hyperterm[0];
            if (show) printf("%i:%e: %e < %e ?  %e\n", i, dblk[0], log(cur.value) + cur.exponent, log(summy.value) + summy.exponent, hyperterm[0]);
            summy.monotonicAdd(cur);
        }
          /*      if ((K > dblk[0] - 0.01)&&(N > dblk[0] - 0.01)) {
            hyperterm[0] *= dblk[0];
            hyperterm[0] *= dblk[0] + N - K - n;
            dblk[0]-= 1.0;
            hyperterm[0] /= -dblk[0] + K;
            hyperterm[0] /= -dblk[0] + n;
            cumNormal(cur, (dblk[0] - dev) / sigma);
            cur *= hyperterm[0];
            if (show) printf("%i: %e < %e ?  %e\n", i, log(cur.value) + cur.exponent, log(summy.value) + summy.exponent, hyperterm[0]);
            summy.monotonicAdd(cur);
        }*/
        if ((dblk[1] + 0.99 < K )&&(dblk[1] + 0.99 < N )){
            hyperterm[1] *= -dblk[1] + K;
            hyperterm[1] *= -dblk[1] + n;
            dblk[1]+= 1.0;
            hyperterm[1] /= dblk[1];
            hyperterm[1] /= dblk[1] + N - K - n;
            cumNormal(cur, (dblk[1] - dev) / sigma);
            cur *= hyperterm[1];
            if (show) printf("%i:%e: %e < %e ?  %e\n", i, dblk[1], log(cur.value) + cur.exponent, log(summy.value) + summy.exponent, hyperterm[0]);
            summy.monotonicAdd(cur);
        }
        old[(i++) & 1] = log(summy.value) + summy.exponent;
        if (i > (N >> 1)) {
            printf("diff = %e, %e + e^%e\n", old[i & 1] - old[(i & 1) ^ 1], summy.value , summy.exponent );
        }
        if (i > N) {
            printf("warning: failed to finish %e %e %i %i %i\n", dev,sigma, n,K,N);
            break;
        }
    }while(old[0] != old[1]);
    summy.exponent += lngammachoose(modeseek, K) + lngammachoose(n - modeseek, N - K) -  lngammachoose(n, N);

    double probsum= -log(summy.value) -summy.exponent;
    if (show) printf("final: %e (%e)\n", (flip) ? log_e_x_m1(probsum): -log_e_x_m1(probsum),-probsum );


     probsum = (flip) ? log_e_x_m1(probsum): -log_e_x_m1(probsum);
    if (ExOp::isValid(probsum)) return(probsum);
    cumNormal(summy, (modeseek - dev) / sigma);
    printf("got NAN... main term predicted is %e + e^%e\n",summy.value, summy.exponent );
return probsum;}
//Series[InverseErfc[2/(1 + E^(1/x))], {x, 0, 2}]
myHashmap<uint32_t, double> getHypergeomericMasses(uint32_t n, uint32_t K, uint32_t N, double probfraction){ myHashmap<uint32_t, double> fout;
    // find Mode
    //if (K == N) {fout[n] = 1.0; return fout;}
    //else if (n == N) {fout[K] = 1.0; return fout;}
    double tot =0;
    if (N < 65536){
        uint32_t mode_p = ((n + 1) *(K +1)) / (N +2);
        double prob_p = exp(lngammachoose(mode_p, K) + lngammachoose(n - mode_p, N - K) -  lngammachoose(n, N));
        tot = prob_p;
        fout[mode_p] = prob_p;
        uint32_t mode_n = mode_p;
        double prob_n = prob_p;
        while(tot < probfraction){
            if (prob_n < prob_p){
                prob_p *= n - mode_p;
                prob_p *= K - mode_p;
                mode_p++;
                prob_p /= mode_p;
                prob_p /= N - K - n + mode_p;
                fout[mode_p] = prob_p;
                tot += prob_p;
            }else{
                if (prob_n == 0.0) { // both extremes must have 0 weight... dont waste your time
                    fout.erase(mode_n);
                    fout.erase(mode_p);
                    return fout;
                }
                prob_n *= mode_n;
                prob_n *= N - K - n + mode_n;
                mode_n--;
                prob_n /= K - mode_n;
                prob_n /= n - mode_n;
                fout[mode_n] = prob_n;
                tot += prob_n;
            }
        }
        if (prob_n == 0.0) fout.erase(mode_n);
        else if (prob_p == 0.0) fout.erase(mode_p);
    }else printf("N is too large.. to do!\n");
return fout;}
void getHypergeomericCover(Tuple<uint32_t> &fout,uint32_t n, uint32_t K, uint32_t N){


}

double HypergeomericLogitPvalBound(uint32_t k, uint32_t n, uint32_t K, uint32_t N){
    uint32_t mode_p = ((n + 1) *(K +1)) / (N +2);
    bool flip;
    if ((flip = (k > mode_p))){
          k = n - k;
          K = N - k;
    }
    double logprob = lngammachoose(k, K) + lngammachoose(n - k, N - K) -  lngammachoose(n, N);
    double firstsum = logprob - log(2);
    double secondsum = log_1_me_x(logprob) - log(2);
return(flip) ? firstsum - secondsum : secondsum - firstsum;}


double HypergeomericNormalPval_slow(double dev, double sigma, uint32_t n, uint32_t K, uint32_t N){
    myHashmap<uint32_t, double> damass = getHypergeomericMasses(n,K,N,0.999999999);
    double sum =0;
    double A = 1 / sigma;
    double B = dev / A;
    for(uint32_t i =0; i < damass.getSize();i++) sum += damass.deref(i) * erf(A*damass.deref_key(i) + B);
return 0.5 + 0.5 * sum;}

/*
SymbolicType::SymbolicType(){}
SymbolicType::~SymbolicType(){}
SymbolicType& SymbolicType::toZero(){dims.toMemfree();return*this;}


ERRCODE SymbolicType::load(FILE *f, unsigned int size){
    ERRCODE fout = ExOp::load(dims,f);
return fout;}
ERRCODE SymbolicType::save(FILE *f) const{
    ERRCODE fout = ExOp::save(dims,f);
return fout;}
*/
//SparseTensor& SparseTensor::toZero(){return *this;}
void SparseTensor::show(FILE*f, int level)const{
    fprintf(f,"SparseTensor\n");
}
uint32_t SymbolicScope::assignAlias(){
    uint32_t r;
    do{ExOp::toRand(r);
    }while(typeflags.find(r) != 0xFFFFFFFF);
    typeflags.addEntry(r);
return r;}
uint32_t SymbolicScope::assignAlias(string name){uint32_t r = this->assignAlias();alias.addEntry(name, r);return(r);}
void SymbolicScope::deleteAlias(uint32_t alias){

}
Tuple<uint32_t> SymbolicScope::mkAliases_routine(Vector<uint32_t> &fout){return Tuple<uint32_t>(fout);}

uint32_t SymbolicScope::prepareIterators(uint32_t alias, uint32_t nbthreads) {
    uint32_t fout = assignAlias();
    //this->members["+iterator"]= fout;
    Tuple<uint32_t> ite; ite.setSize(nbthreads);
    for(int i=0;i<nbthreads;i++) ite[i] = assignAlias();
    uint32_t tlen;
    if (typeflags[currentTarget] < 16){
            switch(typeflags[currentTarget] & 15){
            case 0:
                tlen = int_data[currentTarget].getSize();
                for(int i=0;i<nbthreads;i++){
                    typeflags[ite[i]] = 2;
                    int_ptr[ite[i]] = int_data[currentTarget]() + (i * tlen) / nbthreads;
                    members[ite[i]]["len"] = (((i+1) * tlen) / nbthreads) - ((i * tlen) / nbthreads);
                }
            break;
            case 1:
                tlen = double_data[currentTarget].getSize();
                for(int i=0;i<nbthreads;i++){
                    typeflags[ite[i]] = 3;
                    double_ptr[ite[i]] = double_data[currentTarget]() + (i * tlen) / nbthreads;
                    members[ite[i]]["len"] = (((i+1) * tlen) / nbthreads) - ((i * tlen) / nbthreads);
                }
            break;
            case 2: break;
            case 3: break;
            }

    }

    this->setTarget(fout);
    //(*this) = aliases;
return fout;}
// conditionnal? autoregressive models for Spatial
// template<class I>  TTensor<double,1u,I,TUPLE_FLAG_REMOTE_MEMORY> &target,
bool SymbolicScope::first(uint32_t threadAlias, uint32_t threadInfo){
    uint32_t nb = int_data[threadInfo].getSize();

}
bool SymbolicScope::next(uint32_t threadAlias, uint32_t threadInfo){

}

SymbolicScope& SymbolicScope::operator[](const char* str){
    uint32_t ite = alias.find(str);
    currentTarget = (ite == 0xFFFFFFFF) ? this->assignAlias(str) :  alias.deref(ite);
return *this;}
void SymbolicScope::clearData(uint32_t currentTarget){
    if (typeflags[currentTarget] < 16){
        if (typeflags[currentTarget] & 8) index_data.erase(currentTarget);
        if (typeflags[currentTarget] & 4) dim_data.erase(currentTarget);
        switch(typeflags[currentTarget] & 3){
            case 0: int_data.erase(currentTarget); break;
            case 1: double_data.erase(currentTarget); break;
            case 2: int_ptr.erase(currentTarget); break;
            case 3: double_ptr.erase(currentTarget); break;
        }
    }
    typeflags[currentTarget] =0;
}
SymbolicScope& SymbolicScope::operator=(const TMatrix<double> &data){
    this->clearData(currentTarget);
    typeflags[currentTarget] = 5;
    dim_data.addEntry(currentTarget).setSize(2);
    dim_data[currentTarget][0] = data.sizes[0]; dim_data[currentTarget][1] = data.sizes[1];
    double_data[currentTarget].setSize(data.sizes[0]* data.sizes[1]);
    memcpy(double_data[currentTarget].data , data.data(), sizeof(double) * double_data[currentTarget].getSize());
return *this;}
void SymbolicScope::writeTo(TMatrix<double>& fout)const{
    if (typeflags[currentTarget] < 16){
        switch(typeflags[currentTarget]){
            case 5:
                fout.setSizes(dim_data[currentTarget][0],dim_data[currentTarget][1]);
                memcpy(fout.data(), double_data[currentTarget].data,  sizeof(double) *double_data[currentTarget].getSize());
            break;
        }
    }
}
void SymbolicScope::writeTo(Tuple<double>& fout)const{
    if (typeflags[currentTarget] < 16){
        switch(typeflags[currentTarget]){
            case 1: fout = double_data[currentTarget];break;
        }
    }
}
SymbolicScope& SymbolicScope::operator=(const SparseTuple<double> &data){

return *this;}
SymbolicScope& SymbolicScope::operator=(const Tuple<double> &data){
    this->clearData(currentTarget);
    typeflags[currentTarget] = 1;
    double_data[currentTarget].setSize(data.getSize());
    memcpy(double_data[currentTarget].data , data.data, sizeof(double) * data.getSize());
return *this;}
SymbolicScope& SymbolicScope::operator=(const Tuple<uint32_t> &value){

return *this;}

SymbolicScope& SymbolicScope::operator=(Tuple<double> &&value){

return *this;}
SymbolicScope& SymbolicScope::operator=(Tuple<uint32_t> &&value){

return *this;}
SymbolicScope& SymbolicScope::setTarget(uint32_t _target){
    // assign target an readies iterators
    currentTarget = _target;
return *this;}
SymbolicScope& SymbolicScope::operator=(uint32_t alias){ class Task{
        public:
        SymbolicScope& sym;
        myHashmap<uint32_t, uint8_t*> &target;
        uint32_t blocksize;
        Task(SymbolicScope& _sym, myHashmap<uint32_t, uint8_t*> &_target) : sym(_sym), target(_target){}
        void operator()(uint32_t i,uint32_t maxi){
            uint32_t j;
            /*for(;i < maxi;i++){
                target.heap[i].d = sym.data[sym.currentTarget].heap[i].d;
                target.heap[i].k.k = sym.data[sym.currentTarget].heap[i].k.k;
                target.heap[i].k.d = new uint8_t[blocksize];
                memcpy(target.heap[i].k.d, sym.data[sym.currentTarget].heap[i].k.d, blocksize);
            }*/
        }
    };
   /* Task lotask(*this, data[currentTarget]);
    vars[currentTarget] = vars[alias];
    uint32_t nbblock_mag = (vars[alias].block_bitmag >= 16) ? 0 : 16 - vars[alias].block_bitmag;
    lotask.blocksize = vars[currentTarget].getBlockByteSize();

    uint32_t j;
    uint32_t i = data.getSize();
    lotask.target.heap.setSize(i);
    if (i == 0) return *this;
    i = ((i -1) >> nbblock_mag) << nbblock_mag;
    tb.submit(lotask, i, data.getSize());
    for(;i>0;i=j) {
        j = i - (1 << nbblock_mag);
        tb.submit(lotask, j, i );
    }
    tb.joinAll();*/
return *this;}
SymbolicScope& SymbolicScope::operator+=(uint32_t alias){ class Task{
        public:

        };
return *this;}
SymbolicScope& SymbolicScope::operator-=(uint32_t alias){

return *this;}
SymbolicScope& SymbolicScope::operator*=(uint32_t alias){

return *this;}
SymbolicScope& SymbolicScope::operator/=(uint32_t alias){

return *this;}
SymbolicScope& SymbolicScope::cmp_ll_normal(uint32_t aliasX,uint32_t aliasMu,uint32_t aliasSigma){  class Task{
        public:
        SymbolicScope& sym;
        uint32_t threadNfo;
        Task(SymbolicScope& _sym) : sym(_sym){}
        void operator()(uint32_t threadID){
            uint32_t Xalias = sym.members[threadID]["X"]; // assumed to be 1D
            uint32_t Malias = sym.members[threadID]["Mu"]; // assumed to be 1D
            uint32_t Salias = sym.members[threadID]["Sigma"]; // assumed to be trianglix
            uint32_t Lalias = sym.members[threadID]["LL"]; // assumed to be 1 value

            if (sym.first(threadID,threadNfo)) do{
                myHashmap<uint32_t,double> value;
                uint32_t d = sym.int_data[Xalias][0];
                for(uint32_t i=0;i<d;i++) value[i] = sym.double_ptr[Xalias][i];
                d = sym.int_data[Malias][0];
                for(uint32_t i=0;i<d;i++) value[i] -= sym.double_ptr[Malias][i];

            }while (sym.next(threadID,threadNfo));
        }
    };
    Task lotask(*this);
    lotask.threadNfo = this->prepareIterators(0,0);
}
void SymbolicScope::runFunctor(uint32_t alias){
    // execute function calls
    myHashmap<uint32_t> alias_todo;
    alias_todo.addEntry(alias);
    uint32_t i=0;
    while(i < alias_todo.getSize()){
        printf("dependencies of %i\n", alias);
        ExOp::show(fncs[alias_todo.deref(i)].input_aliases);
        for(uint32_t j=0; j<fncs[alias_todo.deref(i)].input_aliases.getSize();j++){
            if (alias_todo.find(fncs[alias_todo.deref(i)].input_aliases[j]) == 0xFFFFFFFF) alias_todo.addEntry(fncs[alias_todo.deref(i)].input_aliases[j]);
        }
    }
    i= alias_todo.getSize()-1;
    do{
        for(uint32_t j=0;j < fncs[alias_todo.deref(i)].thread_input.getSize();j++){
            this->fncs[alias_todo.deref(i)].fnc(fncs[alias_todo.deref(i)].thread_input[j]);
        }
    }while(i-->0);
}
uint32_t SymbolicScope::assignFunctor( std::function<void(uint32_t)> fnc, Tuple<uint32_t> &input, Tuple<uint32_t> &thread_arg){uint32_t fout = assignAlias();
    ParallelTaskNode& newp = this->fncs[fout];
    newp.input_aliases = input;
    newp.thread_input = thread_arg;
    newp.fnc = fnc;
return fout;}

// that's the thing...
void SymbolicScope::optimize(const char* target){


}


ERRCODE SymbolicScope::load(FILE *f, unsigned int size){
    ERRCODE fout = alias.load(f);
    fout |= members.load(f);
    fout |= typeflags.load(f);
    fout |= int_data.load(f);
    fout |= double_data.load(f);
    fout |= index_data.load(f);
return fout;}
ERRCODE SymbolicScope::save(FILE *f) const{
    ERRCODE fout = alias.save(f);
    fout |= members.save(f);
    fout |= typeflags.save(f);
    fout |= int_data.save(f);
    fout |= double_data.save(f);
    fout |= index_data.save(f);
return fout;}

void SymbolicScope::show(FILE*f, int level){
    /*printf("Symbolic Scope:\n");
    if (auto ite = alias.getIterator()){
        printf("\t%s:", ite().c_str());
    }while(ite++);*/
}

Symbolic& Symbolic::operator=(FunctionNode &fnc){
    fnc.eval(alias);
return(*this);}
Symbolic& Symbolic::operator=(FunctionNode &&fnc){
    fnc.eval(alias);
return(*this);}
Symbolic& Symbolic::operator=(const Tuple<double> &val){

return(*this);}
Symbolic& Symbolic::operator=(Tuple<double> &&val){

return(*this);}

FunctionNode Symbolic::operator+(const Symbolic &other) const{FunctionNode fout(symscp);
    fout.setFunction(LFH_FUNCTION_ADDITION);
    fout.argaliases = symscp.mkAliases(alias, other.alias);
return fout;}
FunctionNode Symbolic::operator-(const Symbolic &other) const{FunctionNode fout(symscp);
    fout.setFunction(LFH_FUNCTION_SUBSTRACTION);
    fout.argaliases = symscp.mkAliases(alias, other.alias);
return fout;}
FunctionNode Symbolic::operator*(const Symbolic &other) const{FunctionNode fout(symscp);
    fout.setFunction(LFH_FUNCTION_MULTIPLICATION);
    fout.argaliases = symscp.mkAliases(alias, other.alias);
return fout;}
FunctionNode Symbolic::operator/(const Symbolic &other) const{FunctionNode fout(symscp);
    fout.setFunction(LFH_FUNCTION_DIVISION);
    fout.argaliases = symscp.mkAliases(alias, other.alias);
return fout;}

void Symbolic::show(FILE* f, int level)const{
    fprintf(f, "Symbolic\n");
}
FunctionNode& FunctionNode::setFunction(LFH_FUNCTION_enum what){func = what;return *this;}
void FunctionNode::setFunction(string name, LFH_FUNCTION_enum what, Tuple<uint32_t> &&input_aliases){
    func = what;
    alias = symscp.assignAlias(name);
    argaliases = input_aliases;
}
void FunctionNode::eval(uint32_t target){
    switch(func){
        case LFH_FUNCTION_ADDITION:
            symscp.setTarget(alias) = argaliases[0];
            for(uint32_t i=1; i < argaliases.getSize();i++) symscp += argaliases[i];
        break;
        case LFH_FUNCTION_SUBSTRACTION:
            symscp.setTarget(alias) = argaliases[0];
            for(uint32_t i=1; i < argaliases.getSize();i++) symscp -= argaliases[i];
        break;
        case LFH_FUNCTION_MULTIPLICATION:
            symscp.setTarget(alias) = argaliases[0];
            for(uint32_t i=1; i < argaliases.getSize();i++) symscp *= argaliases[i];
        break;
        case LFH_FUNCTION_DIVISION:
            symscp.setTarget(alias) = argaliases[0];
            for(uint32_t i=1; i < argaliases.getSize();i++) symscp /= argaliases[i];
        break;
        case LFH_FUNCTION_LL_NORMAL:
            symscp.setTarget(alias) = argaliases[0];
            symscp.cmp_ll_normal(argaliases[1],argaliases[2],argaliases[3]);
        break;
        case LFH_FUNCTION_LL_NEG_BINOMIAL:
            symscp.setTarget(alias) = argaliases[0];
            for(uint32_t i=1; i < argaliases.getSize();i++) symscp /= argaliases[i];
        break;
    }
}
void FunctionNode::backward(uint32_t error_alias){
    switch(func){
        case LFH_FUNCTION_ADDITION:
            symscp.setTarget(alias) = argaliases[0];
            for(uint32_t i=1; i < argaliases.getSize();i++) symscp += argaliases[i];
        break;
    }
}
TensorInnerProduct::TensorInnerProduct(SymbolicScope& _sm): sm(_sm), compute_nbthread(_sm.tb.nbthreads){
    using std::placeholders::_1;
    Tuple<uint32_t> input; input.setSize(3).toZero();
    Tuple<uint32_t> args; args.setSize(_sm.tb.nbthreads);
    for(uint32_t i =0; i < _sm.tb.nbthreads;i++) args[i] =i;
    evalalias = _sm.assignFunctor(std::bind(&TensorInnerProduct::Evaluator::operator(), TensorInnerProduct::Evaluator(*this), _1),input , args);
}
void TensorInnerProduct::Evaluator::operator()(uint32_t id){
    double* buf = new double[trg.tensor.dims[1] * trg.tensor.dims[2]];
    if (auto ite = trg.out_imk(id, trg.compute_nbthread)) do{
        if (auto lite = trg.left[ite()].getIterator()) do{


        }while(lite++);
    }while(ite++);
}
void TensorInnerProduct::run(ThreadBase &tb){
    TensorInnerProduct::Evaluator todo(*this);
    compute_nbthread = tb.nbthreads;

    for(uint32_t i=0;i< compute_nbthread;i++) tb.submit(todo, i);
    tb.joinAll();
}
GaussFunctor::GaussFunctor(SymbolicScope& _sm):sm(_sm), compute_nbthread(_sm.tb.nbthreads){
    //evalalias = _sm.assignFunctor(GaussFunctor::Evaluator(*this));
    using std::placeholders::_1;
    Tuple<uint32_t> input; input.setSize(3).toZero();
    Tuple<uint32_t> args; args.setSize(_sm.tb.nbthreads);
    for(uint32_t i =0; i < _sm.tb.nbthreads;i++) args[i] =i;
    evalalias = _sm.assignFunctor(std::bind(&GaussFunctor::Evaluator::operator(), GaussFunctor::Evaluator(*this), _1),input , args);

}
void GaussFunctor::eval(){
    sm.runFunctor(evalalias);
}

void GaussFunctor::operator()(uint32_t alias){
    if (alias == 0){

    }
}
void GaussFunctor::Evaluator::operator()(uint32_t id){
    Tuple<double> dif;
    if (auto ite = trg.ll_imk(id, trg.compute_nbthread)) do{
        dif = trg.data_acs[ite()] -  trg.mu_acs[ite()];
        auto trig = trg.precisqrt_acs[ite()];
        (*ite) = -0.5 * trig.Xformed_inner_product(dif);
        (*ite) += trig.log_determinant() - 0.91893853320467274178032973640562 * dif.getSize();
        if (ite() == 4){
            dif.show();
            printf("%e is out\n", (*ite));
        }
    }while(ite++);
}
void GaussFunctor::Backward::operator()(uint32_t id){
    Tuple<double> dif;
    if (auto ite = trg.ll_imk(id, trg.compute_nbthread)) do{
        dif = trg.data_acs[ite()] -  trg.mu_acs[ite()];
        auto trig = trg.precisqrt_acs[ite()];
        trg.d_precisqrt_acs[ite()] += trig - Trianglix<double>(dif);
        trg.d_mu_acs[ite()] += trig * dif;
    }while(ite++);
}

void GaussFunctor::run(ThreadBase &tb){
    GaussFunctor::Evaluator todo(*this);
    compute_nbthread =tb.nbthreads;
    for(uint32_t i=0;i< compute_nbthread;i++) tb.submit(todo, i);
    tb.joinAll();
}


void IntegerSampler::init(double* _prob, uint32_t length){
    Tuple<double> buffer;
    imag = 1;
    while(imag < length) imag <<=1;
    prob.setSize(imag);
    imag >>=1;
    buffer.setSize(imag);
    int j;
    for(j = 0; ((j<<1)|1) < length; j++){
        buffer[j] = _prob[(j<<1)] + _prob[(j<<1) |1];
        uint32_t tmp = (uint32_t)(pow(2.0, 32.0f) * _prob[(j<<1) |1] / buffer[j]);
        if ((tmp == 0)&&(_prob[(j<<1)] < _prob[(j<<1) |1]))  prob[j+imag] = 0xFFFFFFFF;
        else prob[j+imag] = tmp;
    }
    if (((j<<1)|1) == length) {prob[j+imag] = 0; buffer[j] =_prob[(j<<1)]; j++; }
    for(;j<imag;j++) buffer[j] = 0;
    j = imag >> 1;
    while(j > 0){
        for(int i = 0; i < j; i++){
            double den = buffer[(i<<1)] + buffer[(i<<1) |1];
            uint32_t tmp = (uint32_t)(pow(2.0, 32.0f) * buffer[(i<<1) |1] / den);
            if ((tmp == 0)&&(buffer[(i<<1)] < buffer[(i<<1) |1]))  prob[i+j] = 0xFFFFFFFF;
            else prob[i+j] = tmp;
            buffer[i] = den;
        }
        j >>= 1;
    }
    //for(int j =0;j < prob.getSize();j++){printf("%i: %X -> %e\n",j, prob[j], pow(2.0, -32.0f) *(double)prob[j] );}
}

uint32_t IntegerSampler::sample()const{

    int i =1;
    uint32_t fout= 0;
    uint32_t j = imag;
    while(i < prob.getSize()){
        uint32_t r;
        do{ExOp::toRand(r);}while(r == 0xFFFFFFFF); // garenties 0 and 1 probabilities can be encoded with 0  0xFFFFFFFF
        if (r < prob[i]){
            fout |= j;
            i = (i << 1) | 1;
        }else i <<= 1;
        j >>= 1;
    }
return fout;}

CMPWindowDistibution& CMPWindowDistibution::setRatesAndPrior(const Tuple<double> &prior, const TMatrix<double> &rates, uint32_t maxlength){
    {
        Tuple<uint32_t,2u> chk; rates.getDims(chk);
        if ((chk[0] != prior.getSize())||(chk[1] != prior.getSize())) {printf("illegal size...\n"); exit(1);}
    }
    condfactor.setSize(maxlength);
    double prob =0;
    log_start_prior.setSize(prior.getSize());
    for(int i=0;i<prior.getSize();i++) prob += prior[i];
    prob = log(prob);
    if (!ExOp::isValid(prob)){
        printf("illegal prior for distribution. (negative value or 0 was provided as probability)\n"); exit(1);
    }
    for(int i=0;i<prior.getSize();i++) log_start_prior[i] = log(prior[i]) - prob;
    transition  =rates;
    for(int j=0;j<prior.getSize();j++){
        prob=0;
        for(int i=0;i<prior.getSize();i++){
            if (i != j) prob += transition(j,i);
        }
        transition(j,j) = prob;
    }
    prob =0;
    for(int i=0;i<prior.getSize();i++){
        prob += exp(log_start_prior[i] -transition(i,i));
    }
    condfactor[0] = log(prob);
    cached_coefs.setSizes(prior.getSize(),prior.getSize());
    for(int j=0;j<prior.getSize();j++){
        for(int i=0;i<prior.getSize();i++){
            if (i ==j) continue;
            cached_coefs(i,j) = exp(log_start_prior[i]) * transition(i,j);
            cached_coefs(i,j) *= transition(j,j) * (transition(j,j) - transition(i,i));
            cached_coefs(i,j) /= transition(j,j) * (1.0 - exp(-transition(i,i))) - transition(i,i) * (1.0- exp(-transition(j,j)));
            if ((!ExOp::isValid(cached_coefs(i,j)))||(cached_coefs(i,j) == 0.0)){
                // did overflow, transition(j,j) too close to transition(i,i)
                cached_coefs(i,j) = exp(log_start_prior[i]) * transition(i,i) * transition(j,j) * transition(i,j) / transition(i,i);
                cached_coefs(i,j) /= 1.0 - exp(-transition(i,i)) * (1.0 + transition(i,i));
            }
        }
    }
}


Tuple< SparseTuple<double> > CMPWindowDistibution::mkSamples(int32_t nbsamples) const{Tuple< SparseTuple<double> > fout;
    fout.setSize(nbsamples);
    double r;
    IntegerSampler mainsm;
    Tuple<double> inp; inp.setSize(log_start_prior.getSize());
    Tuple<IntegerSampler > trsmlr; trsmlr.setSize(log_start_prior.getSize());

    for(int i=0;i<log_start_prior.getSize();i++) inp[i] = exp(log_start_prior[i]);
    mainsm.init(inp(), log_start_prior.getSize());
    for(int i=0;i<log_start_prior.getSize();i++){
        for(int j=0;j<log_start_prior.getSize();j++) inp[j] = (i ==j) ? 0.0 : transition(i,j) / transition(i,i);
        trsmlr[i].init(inp(), log_start_prior.getSize());
    }
    for(int i=0;i<nbsamples;i++){
        uint32_t st = mainsm.sample();
        double r = randExponential() / transition(st,st);

        if (r >= 1.0)  fout[i][st] = 1.0f;
        else{
            uint32_t st2 = trsmlr[st].sample();
            fout[i][st] = r;
            fout[i][st2] = 1.0 -r;
        }
    }
return fout;}



CMPWindowDistibution::LearnScope CMPWindowDistibution::mkEMscope()const{ CMPWindowDistibution::LearnScope fout(*this);
    fout.d_log_start_prior.setSize(log_start_prior.getSize());
    fout.d_transition.setSizes(log_start_prior.getSize(),log_start_prior.getSize());
    fout.lenghtdist.setSize(condfactor.getSize()+1);
return fout;}

void CMPWindowDistibution::LearnScope::init(){
    d_log_start_prior.toZero();
    d_transition.toZero();
    lenghtdist.toZero();
    LL =0.0;
}

double CMPWindowDistibution::LearnScope::EMregist(const SparseTuple<double> &obj){double fout;
    bool islast;
    switch(obj.getSize()){
        case 0: printf("Impossible!\n"); exit(1);
        case 1: if (obj.deref(0) != 1.0) {printf("Impossible!\n"); exit(1);};
            d_log_start_prior[obj.deref_key(0)] += 1.0;
            d_transition(obj.deref_key(0),obj.deref_key(0)) += d_log_1_me_nx_dx(target.transition(obj.deref_key(0),obj.deref_key(0)));
            lenghtdist[0]++;
            return target.log_start_prior[obj.deref_key(0)] + log_1_me_nx(target.transition(obj.deref_key(0),obj.deref_key(0)));
        case 2: {
            fout = target.condfactor[0];
            lenghtdist[1]++;
            islast = (target.condfactor.getSize() == 1);
            Magscaled<double> a;
            Magscaled<double> b;
            int i = obj.deref_key(0);
            int j = obj.deref_key(1);
                a.exponent = target.log_start_prior[i] - target.transition(i,i) * obj.deref(0);
                b.exponent = target.log_start_prior[j] - target.transition(j,j) * obj.deref(1);
                if (!islast){
                    // well not for now but...
                    a.exponent += log_1_me_x(-target.transition(j,j) * obj.deref(1));
                    b.exponent += log_1_me_x(-target.transition(i,i) * obj.deref(0));
                }
                a.value =target.cached_coefs(i,j);
                b.value =target.cached_coefs(j,i);
                a+=b;
                fout = target.condfactor[0] + a.exponent + log(a.value);
            }break;
        default:
            printf("path is too long!\n");
            exit(1);
    }
return fout;}
/*
SparseTuple<double> CMPWindowDistibution::mkSample() const{SparseTuple<double> fout;

    double r;
    IntegerSampler mainsm;
    Tuple<double> inp; inp.setSize(log_start_prior.getSize());
    Tuple<IntegerSampler > trsmlr; trsmlr.setSize(log_start_prior.getSize());

    for(int i=0;i<log_start_prior.getSize();i++) inp[i] = exp(log_start_prior[i]);
    mainsm.init(inp(), log_start_prior.getSize());
    for(int i=0;i<log_start_prior.getSize();i++){
        for(int j=0;j<log_start_prior.getSize();j++) inp[j] = (i ==j) ? 0.0 : transition(i,j) / transition(i,i);
        trsmlr[i].init(inp(), log_start_prior.getSize());
    }


    uint32_t st = mainsm.sample();
    r = randExponential() / transition(st,st);

    if (r >= 1.0) {printf("%i\n",st);fflush(stdout); fout[st] = 1.0f;}
    else{
        uint32_t st2 = trsmlr[st].sample();
        printf("%i -> %i (maybe)\n",st, trsmlr[st].sample());fflush(stdout);
        fout[st] = r;
        fout[st2] = 1.0 -r;
    }

return fout;}*/
double CMPWindowDistibution::operator()(const SparseTuple<double> &target) const{ double fout;
    bool islast;
    switch(target.getSize()){
        case 0: printf("Impossible!\n"); exit(1);
        case 1: if (target.deref(0) != 1.0) {printf("Impossible!\n"); exit(1);};
            return log_start_prior[target.deref_key(0)] + log_1_me_x(-transition(target.deref_key(0),target.deref_key(0)));
        case 2: {
            fout = condfactor[0];
            islast = (condfactor.getSize() == 1);
            Magscaled<double> a;
            Magscaled<double> b;
            int i = target.deref_key(0);
            int j = target.deref_key(1);
                a.exponent = log_start_prior[i] - transition(i,i) * target.deref(0);
                b.exponent = log_start_prior[j] - transition(j,j) * target.deref(1);
                if (!islast){
                    // well not for now but...
                    a.exponent += log_1_me_x(-transition(j,j) * target.deref(1));
                    b.exponent += log_1_me_x(-transition(i,i) * target.deref(0));
                }
                a.value =cached_coefs(i,j);
                b.value =cached_coefs(j,i);
                a+=b;
                fout = condfactor[0] + a.exponent + log(a.value);
            }break;
        default:
            printf("path is too long!\n");
            exit(1);
    }
return fout;}

void CMPWindowDistibution::show(FILE*f, int level)const{
    printf("Continuous Markov Process Window Distribution, with start prior:\n");
    log_start_prior.show();
    printf("Transition rate are:\n");
    transition.show();
    printf("Probability of for vector length is:\n");
    for(int i=0;i< condfactor.getSize();i++) printf("P(L > %i)= %e\n",i+1,exp(condfactor[i]));
}


MPReachDistribution& MPReachDistribution::setDefaultRatesAndPrior(uint32_t nbdims, double pathlength){
    Tuple<double> stprior; stprior.setSize(nbdims).toOne();
    Trianglix<double> trprior; trprior.setSize(nbdims);
    for(int i=0;i< nbdims * (nbdims+1) / 2;i++) trprior.data[i] = 1.0;
return(this->setRatesAndPrior(stprior, trprior, pathlength));}
MPReachDistribution& MPReachDistribution::setRatesAndPrior(const Tuple<double> &prior, const Trianglix<double> &rates, double pathlength){
    if (rates.getSize() != prior.getSize()) {printf("illegal size...\n"); exit(1);}
    if (rates.getSize() == 0) {printf("warning, size is 0...\n"); return *this;}
    double prob =0;
    log_start_prior.setSize(prior.getSize());
    for(int i=0;i<prior.getSize();i++) {
        prob += prior[i];
        if (prior[i] < 0){
            printf("illegal prior for distribution. (negative value  probability)\n"); exit(1);
        }
    }
    prob = log(prob);
    if (!ExOp::isValid(prob)){
        printf("illegal prior for distribution. (negative value  probability)\n"); exit(1);
    }
    for(int i=0;i<prior.getSize();i++) log_start_prior[i] = log(prior[i]) - prob;
    log_transition.setSize(prior.getSize());
    log_transition.data[0] = 0;
    underbound_for_exit_probability = 1.0f - (1.0 / (1.0 + pathlength) );
    int k = 1;
    for(int j=1;j<prior.getSize();j++){
        prob=0;
        for(int i=0;i<j;i++)prob += rates.data[k+i];
        if ((rates.data[k+j] / (prob + rates.data[k+j])) < underbound_for_exit_probability){
            log_transition.data[k+j] = log(underbound_for_exit_probability);
            prob = log(prob) - log(1.0 - underbound_for_exit_probability);
        }else{
            prob += rates.data[k+j];
            prob = log(prob);
            log_transition.data[k+j] = log(rates.data[k+j]) - prob;
        }
        for(int i=0;i<j;i++) log_transition.data[k+i] = log(rates.data[k+i]) - prob;
        k += j+1;
    }

}

double MPReachDistribution::operator()(const SparseTuple<double> &what)const{
    int l = what.getSize();
    if (l == 0) return 0.0;
    Vector<uint32_t> dims; dims.setSize(l);
    for(int i=0;i< l;i++) dims[i] = what.deref_key(i);
    dims.sort();
    l--;
    double LL = log_start_prior[dims[l]];
    for(;l>0;l--){
        LL += log_transition(dims[l], dims[l-1]);
    }
    LL += log_transition(dims[0], dims[0]);
return LL;}
Tuple<SparseTuple<double> > MPReachDistribution::mkSamples(int32_t nbsamples) const{Tuple< SparseTuple<double> > fout;
    fout.setSize(nbsamples);
    if (log_start_prior.getSize() == 0) {printf("warning, not initialized distribution sampled!\n"); return fout;}
    double r;
    Tuple<double> inp; inp.setSize(log_start_prior.getSize());
    Tuple<IntegerSampler > trsmlr; trsmlr.setSize(log_start_prior.getSize());

    for(int i=0;i<log_start_prior.getSize();i++) inp[i] = exp(log_start_prior[i]);
    trsmlr[0].init(inp(), log_start_prior.getSize());

    for(int i=1;i<log_start_prior.getSize();i++){
        inp.setSize(i+1);
        for(int j=0;j<=i;j++) inp[j] = exp(log_transition(i,j));
        trsmlr[i].init(inp(), i+1);
    }
    for(int i=0;i<nbsamples;i++){
        uint32_t st = trsmlr[0].sample();
        uint32_t vfst = st;
        double total =0.0;
        ExOp::toAbs(ExOp::toRand(fout[i][st]));
        total += fout[i][st];
        while(st){
            uint32_t nt = trsmlr[st].sample();
            if (nt == st) break;
            st = nt;
            ExOp::toAbs(ExOp::toRand(fout[i][st]));
            total += fout[i][st];
        }
        if (auto ite = fout[i].getIterator()) do{*ite /= total;} while(ite++);
    }
return fout;}

void MPReachDistribution::show(FILE*f, int level)const{
    printf("Markov Process Reach Distribution, with start log probability prior %i %i:\n", log_start_prior.getSize(), log_transition.getSize() ); fflush(stdout);
    log_start_prior.show();
    printf("Log Transition rate are:\n");
    log_transition.show();
}

Tuple<uint32_t, 2u> MPReachDistribution::findCandidateStep(SparseTuple<double> &fout, const SparseTuple<double> &start_io, uint32_t& outdirhit, uint32_t nbdim, const double *derivative, double LL_prior_scale)const{
    uint32_t nbest = 0xFFFFFFFF;
    Tuple<uint32_t, 2u> suggest;
    suggest[0] = 0xFFFFFFFF;
    double step,tmp;
    double sum =0.0;
    fout.toMemfree();
    int i,l;
    uint32_t ind;



    // Step 1: Constraining to allowed space
    for(i=0;i<nbdim;i++){
        if ((i == 0)||(derivative[i] > derivative[suggest[1]])) suggest[1] = i;
        if ((ind = start_io.find(i)) != 0xFFFFFFFF) {
            fout[i] = derivative[i];
            sum += derivative[i];
            tmp = derivative[i] * start_io.deref(ind);
        }else if ((nbest == 0xFFFFFFFF)||(derivative[i] > derivative[nbest])) nbest = i;
    }

    if ((nbest != 0xFFFFFFFF)&&(derivative[nbest] <= (derivative[nbest] + sum) / (fout.getSize()+1))){ // nothing new to add
        if (fout.getSize() == 1) {// stuck in corner, lets try least terrible direction... who knows?
            fout.deref(0) = -1.0;
            if (nbest == 0xFFFFFFFF) { // should be impossible, unless 1 state in total
                nbest = (rand() % (nbdim - 1));
                if (fout.deref_key(0) <= nbest) nbest++;
            }
            fout[nbest] = 1.0;
            outdirhit = fout.deref_key(0);
            return suggest;
        }
        nbest = 0xFFFFFFFF;
    }

    Vector<uint32_t> dims;
    if (nbest == 0xFFFFFFFF) dims.setSize(fout.getSize());
    else{
        dims.setSize(fout.getSize()+1);
        dims[fout.getSize()] = nbest;
    }
    if (dims.getSize() > 2){
         // try to drop dimension (or the new one if hopeless) ... randomly for now!
        while(true){
            for(int i=0;i< fout.getSize();i++) dims[i] = fout.deref_key(i);
            dims.sort();

            double pivot_val, fact;

            suggest[0] = l = dims.getSize()-1;

            pivot_val = log_start_prior[dims[l-1]] - log_start_prior[dims[l]] - log_transition(dims[l], dims[l-1]);
            for(l-- ; l>0;l--){
                tmp = log_transition(dims[l+1], dims[l-1]) - log_transition(dims[l+1], dims[l]) - log_transition(dims[l], dims[l-1]);
                if (tmp > pivot_val) {suggest[0] = l; pivot_val = tmp;}
            }
            tmp = log_transition(dims[1], dims[1]) - log_transition(dims[1], dims[0]) - log_transition(dims[0], dims[0]);
            if (tmp > pivot_val) {suggest[0] = 0; pivot_val = tmp;}
            suggest[0] = dims[suggest[0]];
            // evaluates if changing sparsity class is worth it
            if (suggest[1] != suggest[0]){
                if (nbest == suggest[0]){
                    if (auto ite = start_io.getIterator()) {
                        tmp = (*ite) * derivative[ite()];
                        while(ite++) tmp += *ite * derivative[ite()];
                    }
                    if ((derivative[nbest] - tmp) <= pivot_val * LL_prior_scale){
                        nbest = 0xFFFFFFFF; // new one would get dropped? do not add it and find a better thing to drop!
                        if (fout.getSize() == 2) {// propose backward motion (as alternative)
                            if (derivative[fout.deref_key(0)] > derivative[fout.deref_key(1)]) {suggest[0] = fout.deref_key(0); suggest[1] = fout.deref_key(1);
                            }else {suggest[0] = fout.deref_key(1); suggest[1] = fout.deref_key(0);}
                            break;
                        }
                        dims.pop_back();
                        continue;
                    }
                    suggest[0] = 0xFFFFFFFF;
                }else if (((derivative[suggest[1]] - derivative[suggest[0]]) * start_io[suggest[0]]) + pivot_val * LL_prior_scale < 0) suggest[0] = 0xFFFFFFFF;
            } else suggest[0] = 0xFFFFFFFF;
            break;
        }
    } else if (fout.getSize() == 2){// propose backward motion (as alternative)
        if (derivative[fout.deref_key(0)] > derivative[fout.deref_key(1)]) {suggest[0] = fout.deref_key(0); suggest[1] = fout.deref_key(1);
        }else {suggest[0] = fout.deref_key(1); suggest[1] = fout.deref_key(0);}
    }

    if (nbest != 0xFFFFFFFF) {sum += derivative[nbest]; fout[nbest] = derivative[nbest];}


    //Step 2: Constraining to unity plane and scale to get maximum step allowed
    sum /= fout.getSize();
    ExOp::toMax(step);
    if (auto ite = fout.getIterator()) do{
        *ite -= sum;
        if ((*ite < 0)&&(step > -start_io[ite()] / (*ite) )) {step = -start_io[ite()] / (*ite); outdirhit = ite();}
    } while(ite++);

    if (auto ite = fout.getIterator()) do{*ite *= step;} while(ite++);
return suggest;}


MPReachDistribution::LearnScope::LearnScope(const MPReachDistribution& _target):target(_target){
    start_count.setSize(_target.log_start_prior.getSize());
    transition_count.setSize(_target.log_start_prior.getSize());
}

void MPReachDistribution::LearnScope::init(){
    start_count.toZero();
    transition_count.toZero();
}
double MPReachDistribution::LearnScope::EMregist(const SparseTuple<double> &what){
    int l = what.getSize();
    Vector<uint32_t> dims; dims.setSize(l);
    for(int i=0;i< l;i++) dims[i] = what.deref_key(i);
    dims.sort();
    l--;
    double LL = target.log_start_prior[dims[l]];
    start_count[dims[l]]++;
    for(;l>0;l--){
        LL += target.log_transition(dims[l], dims[l-1]);
        transition_count(dims[l], dims[l-1])++;
    }
    LL += target.log_transition(dims[0], dims[0]);
    transition_count(dims[0], dims[0])++;
return LL;}

double MPReachDistribution::learn(MPReachDistribution::LearnScope& scp, double prior_count, bool do_show){ double LL=0.0;
    uint32_t tc=0;
    if (do_show){
        printf("stats obtained...\n");
        scp.start_count.show();
        printf("and\n");
        scp.transition_count.show();
    }
    for(int i=0;i < log_start_prior.getSize();i++) tc += scp.start_count[i];
    double lf = log(prior_count + tc);
    for(int i=0;i < log_start_prior.getSize();i++) {
        log_start_prior[i] =  log((prior_count / log_start_prior.getSize()) + scp.start_count[i]) - lf;
        if (scp.start_count[i]) LL += log_start_prior[i] * scp.start_count[i];
    }
    int k =1;
    for(int j=1; j< log_start_prior.getSize(); j++){
        tc =0;
        for(int i=0;i < j;i++) tc += scp.transition_count.data[k+i];
        if (((prior_count / j) + scp.transition_count.data[k+j]) < underbound_for_exit_probability * (prior_count + (scp.transition_count.data[k+j] + tc))){
            // did hit under bound... use the best allowed transitions
            log_transition.data[k+j] = log(underbound_for_exit_probability);
            lf = log((prior_count / j) * (j-1) + tc) - log(1.0 - underbound_for_exit_probability);
        }else{
            lf = log(prior_count + (tc + scp.transition_count.data[k+j]));
            log_transition.data[k+j] = log((prior_count / j) + scp.transition_count.data[k+j]) - lf;
        }
        for(int i=0;i < j;i++) {
            log_transition.data[k+i] = log((prior_count / j) + scp.transition_count.data[k+i]) - lf;
            if (scp.transition_count.data[k+i]) LL += log_transition.data[k+i] * scp.transition_count.data[k+i];
        }
        if (scp.transition_count.data[k+j])  LL += log_transition.data[k+j] * scp.transition_count.data[k+j];
        k+= j+1;
    }

    if (do_show){
        printf("learned log prior:\n");
        log_start_prior.show();
        printf("learned log prior transition:\n");
        log_transition.show();
    }
return LL;}

//double HypergeomericLogitPval(uint32_t k, uint32_t n, uint32_t K, uint32_t N){ // gets the probability masses that accounts for >= probfraction (or all if no threshold provided)


/*
myHashmap<uint32_t, double> getLogErfsum(uint32_t n, uint32_t K, uint32_t N, double probfraction){ myHashmap<uint32_t, double> fout;
    // find Mode
    //if (K == N) {fout[n] = 1.0; return fout;}
    //else if (n == N) {fout[K] = 1.0; return fout;}
    double tot =0;
    if (N < 65536){
        uint32_t mode_p = ((n + 1) *(K +1)) / (N +2);
        double prob_p = exp(lngammachoose(mode_p, K) + lngammachoose(n - mode_p, N - K) -  lngammachoose(n, N));
        tot = prob_p;
        fout[mode_p] = prob_p;
        uint32_t mode_n = mode_p;
        double prob_n = prob_p;
        while(tot < probfraction){
            if (prob_n < prob_p){
                prob_p *= n - mode_p;
                prob_p *= K - mode_p;
                mode_p++;
                prob_p /= mode_p;
                prob_p /= N - K - n + mode_p;
                fout[mode_p] = prob_p;
                tot += prob_p;
            }else{
                if (prob_n == 0.0) { // both extremes must have 0 weight... dont waste your time
                    fout.erase(mode_n);
                    fout.erase(mode_p);
                    return fout;
                }
                prob_n *= mode_n;
                prob_n *= N - K - n + mode_n;
                mode_n--;
                prob_n /= K - mode_n;
                prob_n /= n - mode_n;
                fout[mode_n] = prob_n;
                tot += prob_n;
            }
        }
        if (prob_n == 0.0) fout.erase(mode_n);
        else if (prob_p == 0.0) fout.erase(mode_p);
    }else printf("N is too large.. to do!\n");
return fout;}*/



} // end of namespace

 //   #include "Advanced.cpp"
