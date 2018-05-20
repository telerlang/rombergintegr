#include "rombergintegr.h"
#include <math.h>
#include <stdio.h>

//求T数据序列的公式的累加部分。
//只要算式中有double类型的，就不用把整数输入成小数形式（如2输入成2.0）。因为会自动把整数提升为double类型。
static inline double sumaddend(const struct integration *integr, const int i)
{
	double sum = 0.0L;
	double temp = pow(2.0, i);
	//j从1一直到2的i-1次方。
	for (long j = 1; j <= (long)temp / 2; j++)
		sum += (integr->b - integr->a) / temp * integr->f(integr->a + (2 * j - 1) * (integr->b - integr->a) / temp);
	return sum;
}

//求第1个R，计算出T1,T2,T4,T8,S1,S2,S4,C1,C2,R1。
static void calcfirstR(const struct integration *integr, struct rombergset *rset)
{
	for (int i = 0; i <= 3; i++)
	{
		//计算T0，或T1、T2、T4
		if (i == 0)
			rset->T[i] = (integr->b - integr->a) / 2 * (integr->f(integr->b) + integr->f(integr->a));
		else
			rset->T[i] = rset->T[i - 1] / 2 + sumaddend(integr, i);

		switch (i)
		{
		case 0:
			break;
		case 1: //计算S1--S[1]
			rset->S[i] = (4 * rset->T[i] - rset->T[i - 1]) / 3;
			break;
		case 2: //计算S2--S[2],C1--C[2]
			rset->S[i] = (4 * rset->T[i] - rset->T[i - 1]) / 3;
			rset->C[i] = (16 * rset->S[i] - rset->S[i - 1]) / 15;
			break;
		case 3: //计算S4--S[3],C2--C[3],R1--R[3], delta[3]=R1-C2=R[3]-C[3]
			rset->S[i] = (4 * rset->T[i] - rset->T[i - 1]) / 3;
			rset->C[i] = (16 * rset->S[i] - rset->S[i - 1]) / 15;
			rset->R[i] = (64 * rset->C[i] - rset->C[i - 1]) / 63;
			rset->delta[i] = rset->R[i] - rset->C[i];
			break;
		}
	}
	rset->Tcount = 4; //4个T有效，T0--T[0],T1--T[1],T2--T[2],T4--T[3]
}

//计算Tcount+1的数据序列，并将计算结果填入数组的[Tcount]处。
static inline void calcnext(const struct integration *integr, struct rombergset *rset)
{
	//当前有Tcount个T有效，则T[Tcount-1]是最后一个有效数字，T[Tcount]无效。
	//所以计算下一行是填入T[Tcount]位置。
	int i = rset->Tcount;
	if (i >= RSETMAX)
		return; //要加上这一行，否则会产生内存溢出。
	rset->T[i] = rset->T[i - 1] / 2 + sumaddend(integr, i);
	rset->S[i] = (4 * rset->T[i] - rset->T[i - 1]) / 3;
	rset->C[i] = (16 * rset->S[i] - rset->S[i - 1]) / 15;
	rset->R[i] = (64 * rset->C[i] - rset->C[i - 1]) / 63;
	rset->delta[i] = rset->R[i] - rset->R[i - 1]; //从delta[3]开始数据有效。
	rset->Tcount++;
}

void initrset(struct rombergset *rset)
{
	rset->Tcount = 0;
}

void initintegr(struct integration *integr, double a, double b, double (*f)(double))
{
	integr->a = a;
	integr->b = b;
	integr->f = f;
}

//按T的个数来计算龙贝格数据序列。
//如果输入值大于预设序列数量最大值，或者输入的值小于等于序列的当前数量，则没有计算返回0。
//否则计算序列直到T的数量等于输入count，并返回Tcount值。
int calcrsetbycount(const struct integration *integer, struct rombergset *rset, int count)
{
	if (rset->Tcount == count && count == 0)
	{
		calcfirstR(integer, rset);
		return rset->Tcount;
	}
	if (count > RSETMAX || rset->Tcount >= count)
		return 0;
	if (rset->Tcount < 4)
		calcfirstR(integer, rset);
	while (rset->Tcount < count)
		calcnext(integer, rset);
	return rset->Tcount;
}

//按最小误差来计算龙贝格数据序列
//输出0表示计算到T的最大个数也不能满足最小误差要求。
//输出其它数字表示当前序列的最后一个值（Tcount-1位置）满足最小误差要求。
int calcrsetbydelta(const struct integration *integer, struct rombergset *rset, double mindelta)
{
	if (rset->Tcount < 4)
		calcfirstR(integer, rset);

	int i = 3;
	while (fabs(rset->delta[i]) > mindelta)
	{
		if (rset->Tcount == RSETMAX)
			return 0;
		if (i + 1 == rset->Tcount)
			calcnext(integer, rset);
		i++;
	}
	return i;
}

void printrset(struct rombergset *rset)
{
	if (rset->Tcount == 0)
	{
		printf("is empty set\n");
		return;
	}
	printf("k  powof2    T\t\t S\t\t C\t\t R\t\t delta\n");
	for (int k = 0; k <= rset->Tcount - 1; k++)
	{
		printf("%2d %ld\t ", k, (long)pow(2.0, k));
		printf("%.8f\t %.8f\t %.8f\t %.8f\t %.8f\n", rset->T[k], rset->S[k + 1], rset->C[k + 2], rset->R[k + 3], rset->delta[k + 3]);
	}
}
