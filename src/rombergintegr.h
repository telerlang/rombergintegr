#ifndef _ROMBERGINTEGR_H
#define _ROMBERGINTEGR_H

//求积分函数的上下限和函数体
struct integration
{
	double a;
	double b;
	double (*f)(double);
};

#define RSETMAX 64 //因为长整数类型64位，所以最长定为64应该足够了。

//Tcount存的是T的个数；
//T0存在T[0],T1存在T[1],T4存在T[2];
//S1存在S[1],S2存在S[2],S4存在S[3],S[0]无效；
//C1存在C[2],C2存在C[3],C4存在C[4],C[0]、C[1]无效；
//R1存在R[3],R2存在R[4],R4存在R[5],R[0]、R[1]、R[2]无效；
//delta[4]=R2-R1=R[4]-R[3],delta[5]=R4-R2=R[5]-R[4],delta[0]、delta[1]、delta[2]、delta[3]无效。
//T, S, C, R, delta数据序列的[Tcount]及其之后数据无效。
struct rombergset
{
	int Tcount;
	double T[RSETMAX];
	double S[RSETMAX];
	double C[RSETMAX];
	double R[RSETMAX];
	double delta[RSETMAX];
};

void initrset(struct rombergset *rset);
void initintegr(struct integration *integr, double a, double b, double (*f)(double));
void printrset(struct rombergset *rset);

//按T的个数来计算龙贝格数据序列
//如果输入值大于预设序列数量最大值，或者输入的值小于等于序列的当前数量，则没有计算返回0。
//否则计算序列直到T的数量等于输入count，并返回Tcount值。
int calcrsetbycount(const struct integration *integr, struct rombergset *rset, int Tcountin);

//按最小误差来计算龙贝格数据序列
//输出0表示计算到最大数量也不能满足最小误差要求。
//输出其它数字表示当前序列的最后一个值（Tcount-1位置）满足最小误差要求。
int calcrsetbydelta(const struct integration *integr, struct rombergset *rset, double mindelta);

#endif