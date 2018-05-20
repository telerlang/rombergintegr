#include "rombergintegr.h"
#include <stdio.h>
#include <math.h>

static double f(double x)
{
    if (x == 0)
        return 1;
    else
        return sin(x) / x;
}

void clearset(struct rombergset *rset)
{
    rset->Tcount = 0;
    for (int i = 0; i <= RSETMAX; i++)
    {
        rset->T[i] = 0.0L;
        rset->S[i] = 0.0L;
        rset->C[i] = 0.0L;
        rset->R[i] = 0.0L;
        rset->delta[i] = 0.0L;
    }
}

void testcalcbycount(struct integration *integr, struct rombergset *rset)
{
    int temp;
    //测试输入从0到6，对空序列的计算。
    for (int i = 0; i < 7; i++)
    {
        printf("\ntest by count %d\n", i);
        clearset(rset);
        temp = calcrsetbycount(integr, rset, i);
        printf("return 4(where i<4) or i OK calc to %d\n", temp);
        printrset(rset);
    }

    //当前序列为6，测试输入10。
    printf("\ntest from 6 skip 10\n");
    printf("currect Tcount %d\t", rset->Tcount);
    temp = calcrsetbycount(integr, rset, 10);
    printf("return 10 is OK %d\n", temp);
    printrset(rset);

    //当前序列为10，测试输入8。
    printf("\ncurrect 10, test input 8\n");
    temp = calcrsetbycount(integr, rset, 10);
    printf("return 0 is OK %d\n", temp);
    printf("currect Tcount %d\n", rset->Tcount);

    //测试输入大于RSETMAX的值。
    printf("\ncurrect 10, test input 123456\n");
    temp = calcrsetbycount(integr, rset, 123456);
    printf("return 0 is OK %d\n", temp);
    printf("currect Tcount %d\n", rset->Tcount);
}

void testcalcbydelta(struct integration *integr, struct rombergset *rset)
{
    int temp;

    //测试对1个空序列，输入误差后的计算。
    printf("\ntest empty rset calc by mindelta=0.1e-5\n");
    clearset(rset);
    temp = calcrsetbydelta(integr, rset, 0.1e-5);
    printf("return=Tcount is OK %d \t currect Tcount %d\t delta[Tcount-1] %.22lf\n", temp, rset->Tcount, rset->delta[temp]);
    printrset(rset);

    //测试对序列个数为6个T的，查最小误差为0.1e-1，应该不新计算序列，返回4。
    printf("\ntest 6T rset calc by mindelta=0.1e-1\n");
    clearset(rset);
    calcrsetbycount(integr, rset, 6);
    temp = calcrsetbydelta(integr, rset, 0.1e-1);
    printf("return 4 is OK %d \t currect Tcount %d\t delta[temp] %.22lf\n", temp, rset->Tcount, rset->delta[temp]);
    printrset(rset);

    //测试对序列个数为6个T的，查最小误差为0.1e-14，应该再往后面计算。
    printf("\ntest 6T rset calc by mindelta=0.1e-14\n");
    clearset(rset);
    calcrsetbycount(integr, rset, 6);
    temp = calcrsetbydelta(integr, rset, 0.1e-14);
    printf("return %d \t currect Tcount %d\t delta[temp] %.22lf\n", temp, rset->Tcount, rset->delta[temp]);
    printrset(rset);

    //测试对序列个数为6个T的，查最小误差为0.1e-20，应该再往后面计算。
    printf("\ntest 6T rset calc by mindelta=0.1e-20\n");
    clearset(rset);
    calcrsetbycount(integr, rset, 6);
    temp = calcrsetbydelta(integr, rset, 0.1e-20);
    printf("return %d \t currect Tcount %d\t delta[temp] %.22lf\n", temp, rset->Tcount, rset->delta[temp]);
    printrset(rset);
}

int main()
{
    struct integration integr;
    initintegr(&integr, 0, 1, f);

    struct rombergset rset;
    testcalcbycount(&integr, &rset);
    testcalcbydelta(&integr, &rset);
}