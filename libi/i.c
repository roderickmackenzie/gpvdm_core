// 
// General-purpose Photovoltaic Device Model gpvdm.com- a drift diffusion
// base/Shockley-Read-Hall model for 1st, 2nd and 3rd generation solarcells.
// The model can simulate OLEDs, Perovskite cells, and OFETs.
// 
// Copyright (C) 2012-2017 Roderick C. I. MacKenzie info at gpvdm dot com
// 
// https://www.gpvdm.com
// 
// 
// This program is free software; you can redistribute it and/or modify it
// under the terms and conditions of the GNU Lesser General Public License,
// version 2.1, as published by the Free Software Foundation.
// 
// This program is distributed in the hope it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
// 
// 




/** @file i.c
	@brief Simple functions to read in scientific data from text files and perform simple maths on the data.
*/
#define _FILE_OFFSET_BITS 64
#define _LARGEFILE_SOURCE
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <sim_struct.h>

#include "i.h"
#include "util.h"
#include "cal_path.h"
#include "const.h"
#include <log.h>
#include "inp.h"


static int unused __attribute__((unused));
static char* unused_pchar __attribute__((unused));

void inter_sin(struct istruct *in,long double mag,long double fx,long double delta)
{
int i;
	for (i=0;i<in->len;i++)
	{
		in->data[i]=fabs(mag)*sin(2*PI*fx*(in->x[i]+delta));
	}
}


void inter_add_to_hist(struct istruct* in,long double pos,long double value)
{
int ii=0;
long double min=in->x[0];
long double max=in->x[in->len-1];
long double dx=(max-min)/((long double)in->len);

ii=(int)((pos-min)/dx);

if (ii<in->len)
{
	if (ii>=0)
	{
		in->data[ii]+=value;
	}
}

}

long double inter_get_center_of_peak(struct istruct* in,int i,int window)
{
int delta=window/2;

if ((i-delta<0)||(i+delta>in->len-1))
{
return in->x[i];
}
long double top=0.0;
long double btm=0.0;
int ii;
for (ii=i-delta;ii<i+delta;ii++)
{
top+=in->x[ii]*in->data[ii];
btm+=in->data[ii];
}
return top/btm;
}

/**Find the peaks
@param real real array
@param imag imag array
@param in input data
@param fx frequency
*/
void inter_find_peaks(struct istruct* out,struct istruct* in,int find_max)
{
int i=0;
int ii=0;
int window=2;
long double yn=0.0;
long double yc=0.0;
long double yp=0.0;
int min=0;
int max=0;
int grad_l=0;
int grad_r=0;
for (i=0;i<in->len;i++)
{
		for (ii=1;ii<=window;ii++)
		{
			if ((i-ii)>=0 && ((i-ii)<in->len)) yn=in->data[i-ii];

			if ((i+ii)>=0 &&((i+ii)<in->len)) yp=in->data[i+ii];

			yc=in->data[i];

			if ((yc-yn)<0.0)
			{
				grad_l= -1;
			}else
			{
				grad_l=1;
			}

			if ((yp-yc)<0.0)
			{
				grad_r= -1;
			}else
			{
				grad_r=1;
			}

			if ((grad_l== -1) && (grad_r==1))
			{
				min++;
			}else
			if ((grad_l==1) && (grad_r== -1))
			{
				max++;
			}
		}

		if (min==window)
		{
			if (find_max==FALSE) inter_append(out,inter_get_center_of_peak(in,i,21),in->data[i]);

		}else
		if (max==window)
		{
			if (find_max==TRUE) inter_append(out,inter_get_center_of_peak(in,i,21),in->data[i]);
		}
		min=0;
		max=0;
	}
}

/**Do a DFT
@param real real array
@param imag imag array
@param in input data
@param fx frequency
*/

void inter_dft(long double *real,long double *imag,struct istruct* in,long double fx)
{
long double r=0.0;
long double i=0.0;
int j=0;
long double dt=in->x[1]-in->x[0];
long double len=(long double)in->len;
long double n=(len)*fx*dt;

for (j=0;j<in->len;j++)
{
	r+=in->data[j]*cos(2.0*3.1415926*((long double)j)*n/(len));
	i+=in->data[j]*sin(2.0*3.1415926*((long double)j)*n/(len));
}
*real=r;
*imag=i;
}


int inter_sort_compare(const void *a, const void *b)
{
long double aa = *(long double*)a;
long double bb = *(long double*)b;

if (aa < bb) return -1;
if (aa > bb) return  1;

return 0;
}

/**Do a quick search
@param in input structure
*/
int inter_sort(struct istruct* in)
{
if (in->len==0) return -1;

int i=0;
long double *data=(long double *)malloc(in->len*2*sizeof(long double));

for (i=0;i<in->len;i++)
{
	data[i*2]=in->x[i];
	data[(i*2)+1]=in->data[i];
}

qsort(data, in->len, sizeof(long double)*2, inter_sort_compare);

for (i=0;i<in->len;i++)
{
	in->x[i]=data[i*2];
	in->data[i]=data[(i*2)+1];
}

free(data);
return 0;
}

/**Do a chop search for a value
@param x index array
@param N length
@param find Value to find
*/
int search(long double *x,int N,long double find)
{
if (N==1) return 0;
int pos=N/2;
int step=N/2;
do
{
	step=step/2 + (step % 2 > 0 ? 1 : 0);

	if (x[pos]>find)
	{
		pos-=step;
	}else
	{
		pos+=step;
	}

	if (pos<=0)
	{
		pos=0;
		break;
	}
	if (pos>=(N-1))
	{
		pos=N-1;
		break;
	}
	if (step==0) break;
	if (x[pos]==find) break;
	if ((x[pos]<=find)&&((x[pos+1]>find))) break;

}while(1);

if (pos==(N-1)) pos=N-2;


return pos;
}


/**Get position of quartile
@param in input structure
@param value Value to find
*/
long double inter_get_quartile(struct istruct* in,long double value)
{
int i;
long double sum=0.0;
long double sum2=0.0;
long double dl=0.0;
long double dr=0.0;
long double dx=0.0;
for (i=0;i<in->len;i++)
{
	if (i==0)
	{
		dl=in->x[0];
	}else
	{
		dl=in->x[i-1];
	}

	if (i==in->len-1)
	{
		dr=in->x[in->len-1];
	}else
	{
		dr=in->x[i+1];
	}

	dx=dr-dl;
	sum+=dx*in->data[i];
}


for (i=0;i<in->len;i++)
{
	if (i==0)
	{
		dl=in->x[0];
	}else
	{
		dl=in->x[i-1];
	}

	if (i==in->len-1)
	{
		dr=in->x[in->len-1];
	}else
	{
		dr=in->x[i+1];
	}

	dx=dr-dl;
	sum2+=dx*in->data[i];
	if (sum2>(sum*value)) return in->x[i];
}

return 0.0;
}



/**Translate the input istruct to a log struct
@param in inout istruct
@param out output istruct
*/
void inter_to_log_mesh(struct istruct* out,struct istruct* in)
{
long double a=log10(in->x[0]);
long double b=log10(in->x[in->len-1]);
long double step=(b-a)/((long double)out->len);
int i;
long double pos=a;
for (i=0;i<out->len;i++)
{
out->x[i]=pow(10.0,pos);
out->data[i]=inter_get(in,pow(10.0,pos));
pos+=step;

}

}

/**Use linear interpolation to project an istruct array to a new linear mesh
@param in input istruct
@param out output istruct
*/
void inter_to_new_mesh(struct istruct* in,struct istruct* out)
{
int i;
int ii;
long double pos=in->x[0];
long double delta=(in->x[in->len-1]-in->x[0])/(long double)out->len;
pos+=delta/2.0;
for (i=0;i<out->len;i++)
{
	ii=search(in->x,in->len,pos);

	long double x0=in->x[ii];
	long double x1=in->x[ii+1];

	long double y0=in->data[ii];
	long double y1=in->data[ii+1];

	out->x[i]=pos;
	out->data[i]=y0+((y1-y0)/(x1-x0))*(pos-x0);

	pos+=delta;
}

return;
}


/**Sum a 1D istruct whilst taking the modulus of the data.
@param in input istruct
*/
long double inter_sum_mod(struct istruct* in)
{
int i;
long double sum=0.0;

for (i=0;i<in->len;i++)
{
	sum+=fabs(in->data[i]);
}
return sum;
}

/**Get the average value of the data in a 1D istruct between two points
@param in input istruct
@param start start point
@param stop stop point

*/
long double inter_avg_range(struct istruct* in,long double start,long double stop)
{
int i;
long double sum=0.0;
long double points=0.0;
for (i=0;i<in->len;i++)
{
	if ((in->x[i]>start)&&(in->x[i]<stop))
	{
		sum+=in->data[i];
		points+=1.0;
	}
}
return sum/points;
}

/**Sum a 1D istruct (no modulus)
@param in input istruct
*/
long double inter_sum(struct istruct* in)
{
int i;
long double sum=0.0;

for (i=0;i<in->len;i++)
{
	sum+=in->data[i];
}
return sum;
}

/**Convolve two istructs
@param one input/output istruct
@param two input istruct
*/
void inter_convolve(struct istruct* one,struct istruct* two)
{
int i;
//long double sum=0.0;

for (i=0;i<one->len;i++)
{
	one->data[i]*=two->data[i];
}
}



long double inter_get_diff(char *out_path,struct istruct* one,struct istruct* two,long double start,long double stop,struct istruct* mull)
{
FILE *out;
if (out_path!=NULL) out=fopen(out_path,"w");
long double error=0.0;
if (one->x[0]>start) start=one->x[0];

long double points_max=400.0;
long double points=0.0;
long double dx=(stop-start)/points_max;
long double pos=start;
long double etemp=0.0;
do
{

	if (pos>=start)
	{
		etemp=fabs(inter_get_noend(one,pos)-inter_get_noend(two,pos))*inter_get_noend(mull,pos);
		if (out_path!=NULL) fprintf(out,"%Le %Le\n",pos,etemp);
		error+=etemp;
		points+=1.0;
	}

	pos+=dx;
}while(pos<stop);
if (out_path!=NULL) fclose(out);

return error/points;
}

/**Integrate the data
@param in the structure to integrate
*/
long double inter_intergrate(struct istruct* in)
{
int i;
long double tn=0.0;
long double tp=0.0;
//long double t=0.0;
long double dt=0.0;
//long double Eomega=0.0;
long double sum=0.0;
long double n;

	for (i=0;i<in->len;i++)
	{

		if (i==0)
		{
			tn=in->x[i];
		}else
		{
			tn=in->x[i-1];
		}

		if (i==in->len-1)
		{
			tp=in->x[i];
		}else
		{
			tp=in->x[i+1];
		}

		n=in->data[i];
		dt=fabs((tp-tn)/2.0);

		sum+=n*dt;


	}
return sum;
}

/**Integrate the data between limits
@param in the structure to integrate
@param from lower limit
@param from upper limit
*/
long double inter_intergrate_lim(struct istruct* in,long double from, long double to)
{
int i;
long double tn=0.0;
long double tp=0.0;
//long double t=0.0;
long double dt=0.0;
//long double Eomega=0.0;
long double sum=0.0;
long double n;

	for (i=0;i<in->len;i++)
	{

		if (i==0)
		{
			tn=in->x[i];
		}else
		{
			tn=in->x[i-1];
		}

		if (i==in->len-1)
		{
			tp=in->x[i];
		}else
		{
			tp=in->x[i+1];
		}

		n=in->data[i];
		dt=fabs((tp-tn)/2.0);

		if (tn>from) sum+=n*dt;
		if (tn>to) break;

	}

return sum;
}

long double inter_norm_to_one_range(struct istruct* in,long double start,long double stop)
{
int i;
long double max=0.0;

for (i=0;i<in->len;i++)
{
	if (in->x[i]>start)
	{
		max=in->data[i];
		break;
	}
}

for (i=0;i<in->len;i++)
{
	if ((in->x[i]>start)&&(in->x[i]<stop))
	{
		if (in->data[i]>max) max=in->data[i];
	}
}

for (i=0;i<in->len;i++)
{
in->data[i]/=max;
}

return max;
}

/**Get maximum value of an istruct
@param in input istruct
*/
long double inter_get_max(struct istruct* in)
{
long double max=0.0;

max=inter_get_max_range(in,0,in->len);

return max;
}

long double inter_get_max_range(struct istruct* in,int start, int stop)
{
int i;
long double max=0.0;
if (start<in->len)
{
	max=in->data[start];
}

//if (in->len>0) max=in->data[0];
for (i=start;i<stop;i++)
{
	if (in->data[i]>max) max=in->data[i];
}

return max;
}

int inter_get_max_pos(struct istruct* in)
{
int i;
int pos=0;
long double max=in->data[0];
//if (in->len>0) max=in->data[0];
for (i=0;i<in->len;i++)
{
	if (in->data[i]>max)
	{
		max=in->data[i];
		pos=i;
	}
}

return pos;
}

/**Perform log10 on data in istruct
@param in input istruct
*/
void inter_log_y_m(struct istruct* in)
{
int i;
long double mull=1.0;

for (i=0;i<in->len;i++)
{
	mull=1.0;
	if (in->data[i]<0.0) mull= -1.0;
	in->data[i]=log10(sqrt(in->data[i]*in->data[i]))*mull;
}
}

/**Perform log10 on data in istruct
@param in input istruct
*/
void inter_log_y(struct istruct* in)
{
int i;
for (i=0;i<in->len;i++)
{
	in->data[i]=log10(sqrt(in->data[i]*in->data[i]));
}
}
/**Perform log10 on x axis in istruct
@param in input istruct
*/
void inter_log_x(struct istruct* in)
{
int i;
for (i=0;i<in->len;i++)
{
	if (in->x[i]>0.0) in->x[i]=log10(sqrt(in->x[i]*in->x[i]));
}
}

/**Smooth istruct with a window
@param points input istruct
*/
void inter_smooth_range(struct istruct* out,struct istruct* in,int points,long double x)
{
int i=0;
int ii=0;
int pos=0;
long double tot_point=0.0;
long double tot=0;
	for (i=0;i<in->len;i++)
	{
		for (ii= -points;ii<points+1;ii++)
		{

			pos=i+ii;

			if ((pos<in->len)&&(pos>=0))
			{
				tot+=in->data[pos];//*dx;
				tot_point+=1.0;//dx;
			}
		}

		if (in->x[i]>x)
		{
			out->data[i]=(tot/(long double)tot_point);
		}else
		{
			out->data[i]=in->data[i];
		}
		tot=0.0;
		tot_point=0.0;
	}
}

/**Smooth istruct with a window
@param points input istruct
*/
void inter_smooth(struct istruct* out,struct istruct* in,int points)
{
int i=0;
int ii=0;
int pos=0;
long double tot_point=0.0;
long double tot=0;
	for (i=0;i<in->len;i++)
	{
		for (ii= -points;ii<points+1;ii++)
		{

			pos=i+ii;

			if ((pos<in->len)&&(pos>=0))
			{
				tot+=in->data[pos];//*dx;
				tot_point+=1.0;//dx;
			}
		}

		out->data[i]=(tot/(long double)tot_point);
		tot=0.0;
		tot_point=0.0;
	}
}

/**Remove zeros from the data stored in istruct
@param in input istruct
*/
void inter_purge_zero(struct istruct* in)
{
int i;
int write=0;
int read=0;
for (i=0;i<in->len;i++)
{
	in->data[write]=in->data[read];
	in->x[write]=in->x[read];
	if (in->data[read]==0.0)
	{
		write--;
	}
	read++;
	write++;
}
in->len=write;

inter_realloc(in,in->len);

}

void inter_purge_x_zero(struct istruct* in)
{
int i;
int write=0;
int read=0;
for (i=0;i<in->len;i++)
{
	in->data[write]=in->data[read];
	in->x[write]=in->x[read];

	if (in->x[read]==0.0)
	{
		write--;
	}
	read++;
	write++;
}

in->len=write;

inter_realloc(in,in->len);

}

/**Get the smallest data stored in an istruct array
@param in input istruct
*/
long double inter_get_min(struct istruct* in)
{
int i=0;
long double min=in->data[i];
for (i=0;i<in->len;i++)
{
	if (in->data[i]<min) min=in->data[i];


}
return min;
}

/**Get the smallest data stored in an istruct array
@param in input istruct
*/
long double inter_get_min_range(struct istruct* in,long double min, long double max)
{
int i=0;
long double ret=in->data[i];
for (i=0;i<in->len;i++)
{
	if ((in->x[i]>min)&&(in->x[i]<max))
	{
		if (in->data[i]<ret) ret=in->data[i];
	}

}
return ret;
}

/**Chop an istruct array between two points
@param min min point
@param min max point
*/
void inter_chop(struct istruct* in,long double min, long double max)
{
int i;
int write=0;
int read=0;
for (i=0;i<in->len;i++)
{
	in->data[write]=in->data[read];
	in->x[write]=in->x[read];
	write++;
	if (in->x[read]<min)
	{
		write--;
	}

	if (in->x[read]>max) break;
	read++;




}
in->len=write;

inter_realloc(in,in->len);
}

/**Divide the data in an istruct by a value
@param div value to divide the data by
*/
void inter_div_long_double(struct istruct* in,long double div)
{
int i;
for  (i=0;i<in->len;i++)
{
	in->data[i]/=div;
}

}

/**Rescale the scale and the data
@param in The structure holding the data
@param xmul multiply x axis by this
@param ymul multiply y axis by this
*/

void inter_rescale(struct istruct* in,long double xmul, long double ymul)
{
int i;
for  (i=0;i<in->len;i++)
{
	in->x[i]*=xmul;
	in->data[i]*=ymul;
}

}



/**Add a value from every x element in the array
@param value value to subtract from data
*/
void inter_add_x(struct istruct* in,long double value)
{
int i;
for  (i=0;i<in->len;i++)
{
in->x[i]+=value;
}

}

/**Subtract a value from every data element in the array
@param value value to subtract from data
*/
void inter_sub_long_double(struct istruct* in,long double value)
{
int i;
for  (i=0;i<in->len;i++)
{
in->data[i]-=value;
}

}

/**Divide one array by the other they must be of the same length/x-asis
@param in opperand one, then result
@param in opperand two

*/
void inter_div(struct simulation *sim,struct istruct* one,struct istruct* two)
{
if (one->len!=two->len)
{
	printf_log(sim,"The arrays are not the same length\n");
	exit(0);
}

int i;
for  (i=0;i<one->len;i++)
{
	if (one->x[i]!=two->x[i])
	{
		printf_log(sim,"The arrays do not have the same x axis\n");
		exit(0);
	}
	if (two->data[i]!=0) one->data[i]/=two->data[i];
}

}


/**Add a number to an istruct
@param in input istruct
@param value value to add to istruct
*/
void inter_add_long_double(struct istruct* in,long double value)
{
int i;
for  (i=0;i<in->len;i++)
{
in->data[i]+=value;
}

}

/**Normalize the area under a 1D istruct to one multiplied by a constant
@param in input istruct
@param mul number to multiply the istruct by
*/
void inter_norm_area(struct istruct* in,long double mul)
{
int i;
long double tot=0.0;
long double dx=0.0;
for  (i=0;i<in->len;i++)
{
	if (i==0)
	{
		dx=in->x[1]-in->x[0];
	}else
	if (i==in->len-1)
	{
		dx=in->x[i]-in->x[in->len-2];
	}else
	{
		dx=(in->x[i+1]-in->x[i])/2.0+(in->x[i]-in->x[i-1])/2.0;
	}

	tot+=dx*in->data[i];
}

for  (i=0;i<in->len;i++)
{

	in->data[i]/=tot;
	in->data[i]*=mul;

}

}





int inter_get_col_n(struct simulation *sim,char *name)
{
int i=0;
char temp[10000];
char *token;
int col=0;

FILE *file;
file=fopen(name,"r");
if (file == NULL)
{
	printf_log(sim,"inter_get_col_n can not open file %s\n",name);
	exit(0);
}

do
{
	memset(temp,0,10000);
	unused_pchar=fgets(temp, 10000, file);
	const char s[2] = " ";
	for (i=0;i<strlen(temp);i++)
	{
		if (temp[i]=='\t') temp[i]=' ';
	}

	if ((temp[0]!='#')&&(temp[0]!='\n')&&(temp[0]!='\r')&&(temp[0]!=0))
	{
		col=0;
		token = strtok(temp, s);

		do
		{
			token = strtok(NULL, s);
			if (token==NULL) break;
			if (token[0]!='\n') col++;
		}
		while(token!=NULL);

		col--;
		break;

	}


}while(!feof(file));
fclose(file);
return col;
}

void inter_load_by_col(struct simulation *sim,struct istruct* in,char *name,int col)
{
int i=0;
char temp[1000];
long double x;
long double y;
char *token;
int icol=0;
strcpy(in->name,name);

FILE *file;
file=fopen(name,"r");
if (file == NULL)
{
	printf_log(sim,"inter_load_a can not open file %s\n",name);
	exit(0);
}

inter_init(sim,in);
do
{
	memset(temp,0,1000);
	unused_pchar=fgets(temp, 1000, file);
	const char s[2] = " ";
	for (i=0;i<strlen(temp);i++)
	{
		if (temp[i]=='\t') temp[i]=' ';
	}

	if ((temp[0]!='#')&&(temp[0]!='\n')&&(temp[0]!='\r')&&(temp[0]!=0))
	{
		token = strtok(temp, s);

		sscanf(token,"%Le",&(x));
		if (token != NULL)
		{
			icol=0;
			int ret=0;
			while( token != NULL )
			{
				if (col==icol)
				{
					ret=sscanf(token,"%Le",&(y));
					break;
				}
				token = strtok(NULL, s);
				icol++;
			}

			if (ret==1) inter_append(in,x,y);
		}

	}


}while(!feof(file));
fclose(file);
}


void inter_import_array(struct istruct* in,long double *x,long double *y,int len,int alloc)
{
int i;
in->len=len;

if (alloc==TRUE)
{
inter_alloc(in,in->len);
}

for  (i=0;i<in->len;i++)
{
	in->x[i]=x[i];
	in->data[i]=y[i];
}


}

/**Take the derivative with respect to the x axis of an istruct
@param in input istruct
@param output istruct
*/
void inter_deriv(struct istruct* out,struct istruct* in)
{
int i;
long double yl=0.0;
long double yr=0.0;
long double xl=0.0;
long double xr=0.0;
long double dy=0.0;
for (i=0;i<in->len;i++)
{
if (i==0)
{
	xl=in->x[i];
	yl=in->data[i];
}else
{
	xl=in->x[i-1];
	yl=in->data[i-1];
}

if (i==(in->len-1))
{
	xr=in->x[i];
	yr=in->data[i];
}else
{
	xr=in->x[i+1];
	yr=in->data[i+1];
}
if (yr!=yl)
{
	dy=(yr-yl)/(xr-xl);
}else
{
	dy=0.0;
}
out->x[i]=in->x[i];
out->data[i]=dy;
}

}

/**Invert data on x axis of istruct
@param in istruct to operate on
*/
void inter_swap(struct istruct* in)
{
int i;
long double *xtemp=malloc(sizeof(long double)*in->len);
long double *dtemp=malloc(sizeof(long double)*in->len);


for  (i=0;i<in->len;i++)
{
	dtemp[i]=in->data[i];
	xtemp[i]=in->x[i];
}

for  (i=0;i<in->len;i++)
{
in->data[i]=dtemp[in->len-1-i];
in->x[i]=xtemp[in->len-1-i];
}

free(xtemp);
free(dtemp);
}



int inter_search_token(struct simulation *sim,long double *value,char *token,char *name)
{
int found=FALSE;
char temp[1000];
char temp2[1000];
FILE *file;
file=fopen(name,"r");
if (file == NULL)
{
	printf_log(sim,"inter_load can not open file %s\n",name);
	return -1;
}

do
{
	temp[0]=0;
	gpvdm_fgets(temp, 1000, file);

	if (strcmp_begin(temp,token)==0)
	{
		sscanf(temp,"%s %Le",temp2,value);
		found=TRUE;
		break;
	}


}while(!feof(file));
fclose(file);

if (found==TRUE)
{
	return 0;
}

return -1;
}

/**Load data from a file
@param in the istruct holding the data
@param name The file name.
*/
int inter_load(struct simulation *sim,struct istruct* in,char *name)
{
char temp[1000];
struct inp_file data;
long double x;
long double y;
long len;

strcpy(in->name,name);


inp_init(sim,&data);
if (inp_load(sim,&data,name)!=0)
{
	return -1;
}


inter_init(sim,in);

char *line=inp_get_string(sim,&data);
while(line!=NULL)
{
	//
	//unused_pchar=gpvdm_fgets(temp, 1000, file);
	if ((line[0]!='#')&&(line[0]!='\n')&&(line[0]!='\r')&&(line[0]!=0)&&(line[0]!=0x0D))
	{
		sscanf(line,"%Le %Le",&(x),&(y));

		inter_append(in,x,y);
	}
	line=inp_get_string(sim,&data);
}

inp_free(sim,&data);
return 0;
}

void inter_set_value(struct istruct* in,long double value)
{
int i=0;
for  (i=0;i<in->len;i++)
{
	in->data[i]=value;
}

}

/**Take segments of dx and multiply them by the y-axis.
@param in struct to work on
*/
void inter_y_mul_dx(struct istruct* in)
{
int i=0;

long double dx=0.0;
long double d0=0.0;
long double d1=0.0;
for  (i=0;i<in->len;i++)
{
		if (i==0)
		{
			d0=(in->x[0]);
		}else
		{
			d0=(in->x[i-1]);
		}

		if (i==in->len-1)
		{
			d1=(in->x[i]);
		}else
		{
			d1=(in->x[i+1]);
		}

		dx=(d1-d0)/2.0;
		in->data[i]=in->data[i]*dx;
}

}

/**Make a cumulative graph.
@param in struct to work on
*/
void inter_make_cumulative(struct istruct* in)
{
int i=0;
long double dx=0.0;
long double d0=0.0;
long double d1=0.0;
long double tot=0.0;
for  (i=0;i<in->len;i++)
{
	if (i==0)
	{
		d0=(in->x[0]);
	}else
	{
		d0=(in->x[i-1]);
	}

	if (i==in->len-1)
	{
		d1=(in->x[i]);
	}else
	{
		d1=(in->x[i+1]);
	}

	dx=(d1-d0)/2.0;
	tot+=in->data[i]*dx;
	in->data[i]=tot;
}

}

/**Print istruct to screen
@param in struct to print
*/
void inter_dump(struct simulation *sim,struct istruct* in)
{
int i=0;
for  (i=0;i<in->len;i++)
{
	printf_log(sim,"%Le %Le\n",in->x[i],in->data[i]);
}

}

/**Save an istruct to disk and define path
@param in struct to save
@param path path of output file
@param path name of output file
*/
void inter_save_backup(struct istruct* in,char *name,int backup)
{
char wholename[200];
char backup_file[200];
sprintf(wholename,"%s",name);
if (backup==FALSE)
{
	inter_save(in,wholename);
}else
{

	sprintf(backup_file,"%s.back",name);

	if( access( wholename, F_OK ) != -1 )
	{
	remove(backup_file);
	rename(wholename,backup_file);
	}
	inter_save(in,wholename);
}


}


/**Save an istruct to disk and define path
@param in struct to save
@param path path of output file
@param path name of output file
*/
void inter_save_a(struct istruct* in,char *path,char *name)
{
char wholename[200];
join_path(2, wholename,path,name);
inter_save(in,wholename);
}

void inter_save_seg(struct istruct* in,char *path,char *name,int seg)
{
FILE *file=NULL;
int i=0;

int count_max=in->len/seg;
int count=0;
char temp[1000];
char file_name[1000];
int file_count=0;
for  (i=0;i<in->len;i++)
{
	if (count==0)
	{
		sprintf(file_name,"%s%d.dat",name,file_count);

		join_path(2, temp,path,file_name);

		file=fopen(temp,"w");
		file_count++;
	}
		fprintf(file,"%Le",in->x[i]);
		fprintf(file," %Le",in->data[i]);
	count++;
	fprintf(file,"\n");

	if (count==count_max)
	{
		fclose(file);
		count=0;
	}

}
if (count!=0) fclose(file);

}


/**Save an istruct to disk
@param in struct to save
@param name outputfile
*/
void inter_save(struct istruct* in,char *name)
{
FILE *file;
file=fopen(name,"w");
int i=0;
for  (i=0;i<in->len;i++)
{
	fprintf(file,"%Le %Le\n",in->x[i],in->data[i]);
}

fclose(file);
}

int inter_search_pos(struct istruct* in,long double x)
{
return search(in->x,in->len,x);
}

long double inter_get_raw(long double *x,long double *data,int len,long double pos)
{
long double x0;
long double x1;
long double y0;
long double y1;

long double ret;
int i=0;

if (pos<x[0])
{

return 0.0;
}


if (pos>=x[len-1])
{
	i=len-1;
	x0=x[i-1];
	x1=x[i];

	y0=data[i-1];
	y1=data[i];

}else
{
	i=search(x,len,pos);
	x0=x[i];
	x1=x[i+1];

	y0=data[i];
	y1=data[i+1];
}
ret=y0+((y1-y0)/(x1-x0))*(pos-x0);
return ret;
}

/**Get interpolated data from a data set
@param in The structure holding the data
@param x the position of the data.
@return the interpolated data value
*/
long double inter_get(struct istruct* in,long double x)
{
long double x0;
long double x1;
long double y0;
long double y1;

long double ret;
int i=0;

//if (x>in->x[in->len-1]) return 0.0;
if (x<in->x[0])
{

return 0.0;
}


if (x>=in->x[in->len-1])
{
	i=in->len-1;
	x0=in->x[i-1];
	x1=in->x[i];

	y0=in->data[i-1];
	y1=in->data[i];

}else
{
	i=search(in->x,in->len,x);

	x0=in->x[i];
	x1=in->x[i+1];

	y0=in->data[i];
	y1=in->data[i+1];

}
ret=y0+((y1-y0)/(x1-x0))*(x-x0);
return ret;
}

long double inter_get_hard(struct istruct* in,long double x)
{
//long double x0;
//long double x1;
//long double y0;
//long double y1;

//long double ret;
//int i=0;

if (x>in->x[in->len-1])
{

return 0.0;
}

return inter_get(in,x);
}

long double inter_get_noend(struct istruct* in,long double x)
{

long double x0;
long double x1;
long double y0;
long double y1;

long double ret;
int i=0;

if (x<in->x[0])
{
	return in->data[0];
}

if (x>=in->x[in->len-1])
{
	return in->data[in->len-1];
}


	i=search(in->x,in->len,x);
	x0=in->x[i];
	x1=in->x[i+1];

	y0=in->data[i];
	y1=in->data[i+1];

long double eval=0.0;

if ((y1-y0)==0.0)
{
	eval=0.0;
}else
if ((x-x0)==0.0)
{
	eval=0.0;
}else
{
	eval=((y1-y0)/(x1-x0))*(x-x0);
}

ret=y0+eval;
return ret;
}


void inter_reset(struct istruct* in)
{
in->len=0;
in->max_len=0;
}


long double inter_array_get_max(long double *data,int len)
{
int i;
long double max=data[0];
for (i=0;i<len;i++)
{
if (max<data[i]) max=data[i];
}
return max;
}


void inter_join_bins(struct istruct* in,long double delta)
{
int i;
long double tot=0.0;
int pos=0;
long double bin=in->x[0];
int move_on=FALSE;
for (i=0;i<in->len;i++)
{
	move_on=FALSE;

	if (fabs(bin-in->x[i])<delta)
	{
		tot+=in->data[i];

	}else
	{
		move_on=TRUE;
	}

	if (i==in->len-1)
	{
		move_on=TRUE;
	}

	if (move_on==TRUE)
	{
		in->data[pos]=tot;
		in->x[pos]=bin;
		bin=in->x[i];
		tot=in->data[i];
		pos++;

	}

}
in->len=pos;
return;
}



