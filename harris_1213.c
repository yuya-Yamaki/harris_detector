//harris_1208
//0-255の制限なし
//サイトのプログラムはIx^2,Iy^2,IxIy^2の画像出力にあたりuint8に正規化していたが，
// 計算上は制約をかけておらずそのまま計算していた．
// →よって同様に制約をかけずに同様の正規化と出力を行う．
// →また，ハリスによるヘッセ行列ではガウシアンフィルタをかけたもので定義を行い，固有値を求める

//harris_another_from_github.c
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#define lambda 0.04
#define th 0.0001
#define TOL 1.0e-10
#define N 2
#define MAX 100

typedef struct Image
{
    int width, height;
    int maxval;
    uint8_t **val;
} Image;

typedef struct Image_harris
{
    int width, height;
    int maxval;
    double **val;
} Image_harris;

typedef struct edge_corner
{
  int y,x;
}edge_corner;


/****************************画像の入出力，メモリの確保*****************************/
void alloc_mem(Image *img, int *img_data)
{
    int i;

    img->val = (uint8_t **)calloc(img_data[0], sizeof(uint8_t *));
    for (i = 0; i < img_data[0]; i++)
    {
        img->val[i] = (uint8_t *)calloc(img_data[1], sizeof(uint8_t));
    }
}

void alloc_mem_harris(Image_harris *img, int *img_data)
{
    int i;

    img->val = (double **)calloc(img_data[0], sizeof(double *));
    for (i = 0; i < img_data[0]; i++)
    {
        img->val[i] = (double *)calloc(img_data[1], sizeof(double));
    }
}

void free_image(Image *img)
{
    if (img != NULL && img->val != NULL)
    {
        free(img->val);
        free(img);
    }
    else
    {
        fprintf(stderr, "! error in free_image()\n");
        exit(1);
    }
}

void free_harris_image(Image_harris *img)
{
    if (img != NULL && img->val != NULL)
    {
        free(img->val);
        free(img);
    }
    else
    {
        fprintf(stderr, "! error in free_image()\n");
        exit(1);
    }
}

void read_img(char *filename, Image *org)
{
    int width, height;
    int maxval;
    char tmp[256];
    int i, j;
    int img_data[3];
    FILE *fp = fopen(filename, "rb");
    fgets(tmp, 256, fp);
    if (tmp[0] != 'P' || tmp[1] != '5')
    {
        fprintf(stderr, "Not a PGM file!\n");
        exit(1);
    }

    while (*(fgets(tmp, 256, fp)) == '#')
        ;
    sscanf(tmp, "%d %d", &width, &height);
    while (*(fgets(tmp, 256, fp)) == '#')
        ;
    sscanf(tmp, "%d", &maxval);
    org->width = width;
    org->height = height;
    org->maxval = maxval;
    img_data[0] = height;
    img_data[1] = width;
    img_data[2] = maxval;
    alloc_mem(org, img_data);

    for (i = 0; i < height; i++)
    {
        for (j = 0; j < width; j++)
        {
            org->val[i][j] = (uint8_t)fgetc(fp);
        }
    }
    fclose(fp);
}

Image *normalize_for_write_pgm(double **org, int *img_data)
{
    int i, j;
    double max, min;
    Image *mod = (Image *)calloc(1, sizeof(Image));
    alloc_mem(mod, img_data);
    //動的メモリの確保

    max = org[0][0];
    min = org[0][0];
    for (i = 0; i < img_data[0]; i++)
    {
        for (j = 0; j < img_data[1]; j++)
        {
            if (max < org[i][j])
                max = org[i][j];
            if (min > org[i][j])
                min = org[i][j];
        }
    }
    min = fabs(min);
    //printf("max = %f, min = %f\n", max, min);
    for (i = 0; i < img_data[0]; i++)
    {
        for (j = 0; j < img_data[1]; j++)
        {
            mod->val[i][j] = (uint8_t)(255.0 * ((org[i][j] - min) / (max - min)));
        }
    }
    return mod;
}

void write_pgm(Image *img, char *filename, int *img_data)
{
    int i, j;
    FILE *fp;
    fp = fopen(filename, "wb");
    fprintf(fp, "P5\n%d %d\n%d\n", img_data[1], img_data[0], img_data[2]);
    for (i = 0; i < img_data[0]; i++)
    {
        for (j = 0; j < img_data[1]; j++)
        {
            putc(img->val[i][j], fp);
        }
    }
    fclose(fp);
    return;
}

void write_pgm_harris(Image *org, Image_harris *img, char *filename)
{
    int i, j;
    FILE *fp;
    int tmp = 0;
    fp = fopen(filename, "wb");
    fprintf(fp, "P5\n%d %d\n%d\n", img->width, img->height, img->maxval);
    for (i = 0; i < img->height; i++)
    {
        for (j = 0; j < img->width; j++)
        {
            if (img->val[i][j] == 255)
            {
                tmp = (int)(img->val[i][j]);
                putc(tmp, fp);
            }
            else
            {
                tmp = org->val[i][j];
                putc(tmp, fp);
            }
        }
    }
    fclose(fp);
    return;
}

void write_yuv_harris(Image_harris *img, Image_harris *img_edge, Image *org, char *filename)
{
    int i, j;
    FILE *fp;
    int yuv = 0;
    int tmp = 0;

    fp = fopen(filename, "wb");
    for (yuv = 0; yuv < 3; yuv++)
    {
        for (i = 0; i < img->height; i++)
        {
            for (j = 0; j < img->width; j++)
            {

                switch (yuv)
                {
                case 0:
                    tmp = org->val[i][j];
                    putc(tmp, fp);
                    break;
                case 2:
                    if (img->val[i][j] == 0)
                        tmp = 128;
                    else
                        tmp = (uint8_t)img->val[i][j];
                    putc(tmp, fp);
                    break;
                case 1:
                    if (img_edge->val[i][j] == 0)
                        tmp = 128;
                    else
                        tmp = 255;
                    putc(tmp, fp);
                    break;
                } //switch
            }
        } //走査
    }     //yuv
    fclose(fp);
    return;
}
/********************************************************************************************/


/**********************************************************************************************
 *
 *
 *              フィルタ処理，その他計算，ハリスのコーナー値の算出
 *
 *
************************************************************************************************/
//convolve
double **convolve(Image *org, int *fil, int *img_data)
{
    int x, y;
    int i, j;
    int s, t;
    int z = 0;
    double tmp = 0;
    double **mod;
    //動的メモリの確保
    mod = (double **)calloc(img_data[0], sizeof(double *));
    for (i = 0; i < img_data[0]; i++)
    {
        mod[i] = (double *)calloc(img_data[1], sizeof(double));
    }
    printf("start convolve function\n");
    //走査
    for (i = 0; i < img_data[0]; i++)
    {
        for (j = 0; j < img_data[1]; j++)
        {
            //フィルタの畳み込み処理
            for (s = -1; s <= 1; s++)
            {
                for (t = -1; t <= 1; t++)
                {
                    x = j + t;
                    if (x < 0)
                        x = 0;
                    else if (x > img_data[1] - 1)
                        x = img_data[1] - 1;
                    y = i + s;
                    if (y < 0)
                        y = 0;
                    else if (y > img_data[0] - 1)
                        y = img_data[0] - 1;

                    tmp += org->val[y][x] * fil[z];
                    z++;
                }
            }
            mod[i][j] = tmp;
            z = 0;
            tmp = 0;
        }
    }
    printf("finish convolve function\n");
    return mod;
}
/*
Image *lowpass_org(Image *org, int *img_data)
{
    int x, y;
    int i, j;
    int s, t;
    double tmp = 0;
    Image *mod = (Image *)calloc(1, sizeof(Image));
    alloc_mem(mod, img_data);
    printf("start convolve function\n");
    //走査
    for (i = 0; i < img_data[0]; i++)
    {
        for (j = 0; j < img_data[1]; j++)
        {
            //フィルタの畳み込み処理
            for (s = -1; s <= 1; s++)
            {
                for (t = -1; t <= 1; t++)
                {
                    x = j + t;
                    if (x < 0)
                        x = 0;
                    else if (x > img_data[1] - 1)
                        x = img_data[1] - 1;
                    y = i + s;
                    if (y < 0)
                        y = 0;
                    else if (y > img_data[0] - 1)
                        y = img_data[0] - 1;

                    tmp += org->val[y][x];
                }
            }
            tmp = tmp / 9;
            if (tmp > 255)
                tmp = 255;
            else if (tmp < 0)
                tmp = 0;
            mod->val[i][j] = tmp;
            tmp = 0;
        }
    }
    printf("finish convolve function\n");
    return mod;
}
*/
Image *lowpassGauss_org(Image *org, double *fil, int *img_data)
{
    int x, y;
    int i, j;
    int s, t;
    int z = 0;
    double tmp = 0;
    Image *mod = (Image *)calloc(1, sizeof(Image));
    alloc_mem(mod, img_data);
    printf("start convolve function\n");
    //走査
    for (i = 0; i < img_data[0]; i++)
    {
        for (j = 0; j < img_data[1]; j++)
        {
            //フィルタの畳み込み処理
            for (s = -1; s <= 1; s++)
            {
                for (t = -1; t <= 1; t++)
                {
                    x = j + t;
                    if (x < 0)
                        x = 0;
                    else if (x > img_data[1] - 1)
                        x = img_data[1] - 1;
                    y = i + s;
                    if (y < 0)
                        y = 0;
                    else if (y > img_data[0] - 1)
                        y = img_data[0] - 1;

                    tmp += org->val[y][x] * fil[z];
                    z++;
                }
            }
            mod->val[i][j] = tmp;
            z = 0;
            tmp = 0;
        }
    }
    printf("finish convolve Gauss_lowpass function\n");
    return mod;
}

//convolve_Gaussian
double **convolve_Gaus(double **org, double *fil, int *img_data)
{
    int x, y;
    int i, j;
    int s, t;
    int z = 0;
    double tmp = 0;
    double **mod;
    //動的メモリの確保
    mod = (double **)calloc(img_data[0], sizeof(double *));
    for (i = 0; i < img_data[0]; i++)
    {
        mod[i] = (double *)calloc(img_data[1], sizeof(double));
    }

    for (i = 0; i < img_data[0]; i++)
    {
        for (j = 0; j < img_data[1]; j++)
        {
            //Gaussian畳み込み
            for (s = -1; s <= 1; s++)
            {
                for (t = -1; t <= 1; t++)
                {
                    x = j + t;
                    if (x < 0)
                        x = 0;
                    else if (x > img_data[1] - 1)
                        x = img_data[1] - 1;
                    y = i + s;
                    if (y < 0)
                        y = 0;
                    else if (y > img_data[0] - 1)
                        y = img_data[0] - 1;

                    tmp += org[y][x] * fil[z];
                    z++;
                }
            }

            mod[i][j] = tmp;
            z = 0;
            tmp = 0;
        }
    }
    return mod;
}


double **square_over(double **input, int *img_data)
{
    int i, j;
    double tmp;
    double **mod;
    //動的メモリの確保
    mod = (double **)calloc(img_data[0], sizeof(double *));
    for (i = 0; i < img_data[0]; i++)
    {
        mod[i] = (double *)calloc(img_data[1], sizeof(double));
    }

    for (i = 0; i < img_data[0]; i++)
    {
        for (j = 0; j < img_data[1]; j++)
        {
            tmp = input[i][j] * input[i][j];
            mod[i][j] = tmp;
        }
    }

    return mod;
}

double **dxdy_calc(double **dx, double **dy, int *img_data)
{
    int i, j;
    double **dxdy;
    //動的メモリの確保
    dxdy = (double **)calloc(img_data[0], sizeof(double *));
    for (i = 0; i < img_data[0]; i++)
    {
        dxdy[i] = (double *)calloc(img_data[1], sizeof(double));
    }

    for (i = 0; i < img_data[0]; i++)
    {
        for (j = 0; j < img_data[1]; j++)
        {
            dxdy[i][j] = dx[i][j] * dy[i][j];
        }
    }

    return dxdy;
}

double **dx_plux_dy_calc(double **g_dx2, double **g_dy2, int *img_data)
{
    int i, j;
    double **g_dx2_plus_g_dy2;
    //動的メモリの確保
    g_dx2_plus_g_dy2 = (double **)calloc(img_data[0], sizeof(double *));
    for (i = 0; i < img_data[0]; i++)
    {
        g_dx2_plus_g_dy2[i] = (double *)calloc(img_data[1], sizeof(double));
    }

    for (i = 0; i < img_data[0]; i++)
    {
        for (j = 0; j < img_data[1]; j++)
        {
            g_dx2_plus_g_dy2[i][j] = g_dx2[i][j] + g_dy2[i][j];
        }
    }

    return g_dx2_plus_g_dy2;
}

Image_harris *harris_calc(double **g_dx2, double **g_dy2, double **g_dxdy2, double **g_dx2_plus_g_dy2_sqr, int *img_data)
{
    int i, j;

    Image_harris *harris = (Image_harris *)calloc(1, sizeof(Image_harris));
    alloc_mem_harris(harris, img_data);

    for (i = 0; i < img_data[0]; i++)
    {
        for (j = 0; j < img_data[1]; j++)
        {
            harris->val[i][j] = (g_dx2[i][j] * g_dy2[i][j] - g_dxdy2[i][j]) - lambda * g_dx2_plus_g_dy2_sqr[i][j];
            //harris = (AB - C^2) - k((A+B)^2)
        }
    }

    harris->height = img_data[0];
    harris->width = img_data[1];
    harris->maxval = img_data[2];
    return harris;
}

void harris_feature(Image_harris *harris, Image_harris *harris_edge, int *img_data)
{
    int i, j;
    double max, min;
    min = harris->val[0][0];
    max = harris->val[0][0];
    int NumOfPel = img_data[0] * img_data[1];
    edge_corner *edge, *corner;
    edge = (edge_corner *)calloc(NumOfPel, sizeof(edge_corner));
    corner = (edge_corner *)calloc(NumOfPel, sizeof(edge_corner));
    int count_edge = 0, count_corner = 0;
    for (i = 0; i < img_data[0]; i++)
    {
        for (j = 0; j < img_data[1]; j++)
        {
            if (max < harris->val[i][j])
                max = harris->val[i][j];
            if (min > harris->val[i][j])
                min = harris->val[i][j];
        }
    }
    printf("max = %f, min = %f\n", max, min); //max = 5686533083.800625, min = -4380401328.090000
    for (i = 0; i < img_data[0]; i++)
    {
        for (j = 0; j < img_data[1]; j++)
        {
            if (harris->val[i][j] >= th * max){
                harris->val[i][j] = 255;
                harris_edge->val[i][j] = 0;
                edge[count_edge].y = i;
                edge[count_edge].x = j;
                count_edge++;
            } else if (harris->val[i][j] <= (-1) * th * max){
                harris_edge->val[i][j] = 255;
                harris->val[i][j] = 0;
                corner[count_corner].y = i;
                corner[count_corner].x = j;
                count_corner++;
            } else {
                harris->val[i][j] = 0;
                harris_edge->val[i][j] = 0;
            }
        }
    }
    printf("count_edge = %d, count_corner = %d\n", count_edge, count_corner);
}

/*eigen val by jacobi method*/

/*********  固有値をJacobi法による求める副関数  ************************/
/*   戻り値（エラーコード）：  0=正常,      9: 異常                    */
/*   入力： a[n×n]=行列，  n: 要素数                                  */
/*   出力： 固有値（解）は a の対角要素                                */
/*          固有ベクトルは x[n×n]                                     */
/************************************************************************/
int S_Jacobi(double *a, double *x, int n)
{
    int i, j, k, m, count, status;
    double amax, amax0, theta, co, si, co2, si2, cosi, pi = 4.0 * atan(1.0);
    double aii, aij, ajj, aik, ajk;

    //   初期値設定
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (i == j)
                x[n * i + j] = 1.0;
            else
                x[n * i + j] = 0.0;
        }
    }

    //   反復計算
    count = 0;
    status = 9;
    while (count <= MAX)
    {
        //  非対角要素の最大値を探索
        amax = 0.0;
        for (k = 0; k < n - 1; k++)
        {
            for (m = k + 1; m < n; m++)
            {
                amax0 = fabs(a[n * k + m]);
                if (amax0 > amax)
                {
                    i = k;
                    j = m, amax = amax0;
                }
            }
        }
        //  収束判定
        if (amax <= TOL)
        {
            status = 0;
            break;
        }
        else
        {
            aii = a[n * i + i];
            aij = a[n * i + j];
            ajj = a[n * j + j];
            //   回転角度計算
            if (fabs(aii - ajj) < TOL)
            {
                theta = 0.25 * pi * aij / fabs(aij);
            }
            else
            {
                theta = 0.5 * atan(2.0 * aij / (aii - ajj));
            }
            co = cos(theta);
            si = sin(theta);
            co2 = co * co;
            si2 = si * si;
            cosi = co * si;

            //   相似変換行列
            a[n * i + i] = co2 * aii + 2.0 * cosi * aij + si2 * ajj;
            a[n * j + j] = si2 * aii - 2.0 * cosi * aij + co2 * ajj;
            a[n * i + j] = 0.0;
            a[n * j + i] = 0.0;
            for (k = 0; k < n; k++)
            {
                if (k != i && k != j)
                {
                    aik = a[n * k + i];
                    ajk = a[n * k + j];
                    a[n * k + i] = co * aik + si * ajk;
                    a[n * i + k] = a[n * k + i];
                    a[n * k + j] = -si * aik + co * ajk;
                    a[n * j + k] = a[n * k + j];
                }
            }

            //   固有ベクトル
            for (k = 0; k < n; k++)
            {
                aik = x[n * k + i];
                ajk = x[n * k + j];
                x[n * k + i] = co * aik + si * ajk;
                x[n * k + j] = -si * aik + co * ajk;
            }
            count++;
        }
        // printf(" S_Jacobi> iter=%d", count);
        // for (i = 0; i < n; i++)
        //     printf(" %10.6f,", a[n * i + i]);
        // printf("\n");
    }
    return status;
}
/*********  行列を表示する副関数  *******************************/
/*       入力；   a= 行列, n= 表示する列数, m= 表示する行数     */
/****************************************************************/
void S_OutMat(double *a, int n, int m)
{
    int i, j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < m; j++)
            printf("  %10.6f", a[i * m + j]);
        printf("\n");
    }
}

double **eigenval(double **g_dx2, double **g_dy2, double **g_dxdy, int *img_data)
{
    double **mod;
    double Hessian[4];
    int cnt;
    int i, j, s, status;
    double x[N * N];
    FILE *fp;
    fp = fopen("eigen_val.txt", "wb");
    //動的メモリの確保
    mod = (double **)calloc(img_data[0], sizeof(double *));
    for (i = 0; i < img_data[0]; i++)
    {
        mod[i] = (double *)calloc(img_data[1], sizeof(double));
    }

    for (i = 0; i < img_data[0]; i++)
    {
        for (j = 0; j < img_data[1]; j++)
        {
            Hessian[0] = g_dx2[i][j];
            Hessian[1] = g_dxdy[i][j];
            Hessian[2] = g_dxdy[i][j];
            Hessian[3] = g_dy2[i][j];
            //printf(" Jacobi法による固有値解析 \n");
            status = S_Jacobi(Hessian, x, N);
            //Mean Curvature:H
            //mod[i][j] = (Hessian[0] + Hessian[3]) / 2;
            //Gaussian Curvature:K
            mod[i][j] = (Hessian[0] * Hessian[3]);

            //固有値の確認
            fprintf(fp, "%7.2f,\t%7.2f,\t%7.2f,\t%7.2f\n", Hessian[0], Hessian[1], Hessian[2], Hessian[3]);

            //   解析結果の出力
            // if (status == 0)
            // {
            //     printf(" 固有値 \n");
            //     for (s = 0; s < N; s++)
            //         printf(" %10.6f", Hessian[N * s + s]);
            //     printf("\n");
            //     printf(" 固有ベクトル \n");
            //     S_OutMat(x, N, N);
            // }
        }
    }
    fclose(fp);
    return mod;
}

/**************************************************************************/


/******************************************************************************************************
 *
 *
 *
 *                            main関数
 *
 *
 *
*******************************************************************************************************/
int main(int argc, char *argv[])
{
    int i;
    int img_data[3];
    double **dx, **dy, **dx2, **dy2, **dxdy,
        **g_dx2, **g_dy2, **g_dxdy, **g_dxdy2,
        **g_dx2_plus_g_dy2, **g_dx2_plus_g_dy2_square;

    int sobel_x[9] = {
        -1, 0, 1,
        -2, 0, 2,
        -1, 0, 1};

    int sobel_y[9] = {
        -1, -2, -1,
        0, 0, 0,
        1, 2, 1};

    double Gaussian[9] = {
        0.0625, 0.125, 0.0625,
        0.125, 0.25, 0.125,
        0.0625, 0.125, 0.0625};

    if (argc < 3)
    {
        fprintf(stderr, "Usage: %s input.pgm output.pgm\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    Image *org = (Image *)calloc(1, sizeof(Image));
    read_img(argv[1], org);
    img_data[0] = org->height;
    img_data[1] = org->width;
    img_data[2] = org->maxval;

    Image *org_lowpass = (Image *)calloc(1, sizeof(Image));
    Image_harris *harris = (Image_harris *)calloc(1, sizeof(Image_harris));
    Image_harris *harris_edge = (Image_harris *)calloc(1, sizeof(Image_harris));
    alloc_mem_harris(harris_edge, img_data);

    //原画像に対し処理前に平滑化
    org_lowpass = lowpassGauss_org(org, Gaussian, img_data);
    printf("convolve sobel dx\n");
    dx = convolve(org_lowpass, sobel_x, img_data);
    printf("convolve sobel dy\n");
    dy = convolve(org_lowpass, sobel_y, img_data);

    //dx2, dy2, dxdyの算出
    printf("convolve sobel dx2\n");
    dx2 = square_over(dx, img_data);
    Image *dx2_for_write;
    dx2_for_write = normalize_for_write_pgm(dx2, img_data);
    write_pgm(dx2_for_write, "check_dx2.pgm", img_data);
    printf("convolve sobel dy2\n");
    dy2 = square_over(dy, img_data);
    Image *dy2_for_write;
    dy2_for_write = normalize_for_write_pgm(dy2, img_data);
    write_pgm(dy2_for_write, "check_dy2.pgm", img_data);
    printf("convolve sobel dxdy\n");
    dxdy = dxdy_calc(dx, dy, img_data);
    Image *dxdy_for_write;
    dxdy_for_write = normalize_for_write_pgm(dxdy, img_data);
    write_pgm(dxdy_for_write, "check_dxdy.pgm", img_data);

    //ガウシアンフィルタ処理（二階微分と同義）
    printf("convolve sobel g_dx2\n");
    g_dx2 = convolve_Gaus(dx2, Gaussian, img_data); //A
    printf("convolve sobel g_dy2\n");
    g_dy2 = convolve_Gaus(dy2, Gaussian, img_data); //B
    printf("convolve sobel g_dxdy\n");
    g_dxdy = convolve_Gaus(dxdy, Gaussian, img_data); //C

    //行列Hの固有値の算出（jacobi法）
    double **eigen;
    eigen = eigenval(g_dx2, g_dy2, g_dxdy, img_data);
    // Image *eigen_img;
    // eigen_img = normalize_for_write_pgm(eigen, img_data);
    // write_pgm(eigen_img, "Gauss_curvature.pgm", img_data);

    printf("convolve sobel g_dxdy2\n");
    g_dxdy2 = square_over(g_dxdy, img_data); //C^2
    printf("convolve sobel g_dx2_plus_g_dy2\n");
    g_dx2_plus_g_dy2 = dx_plux_dy_calc(g_dx2, g_dy2, img_data); //A+B
    printf("convolve sobel g_dx2_plus_g_dy2_square\n");
    g_dx2_plus_g_dy2_square = square_over(g_dx2_plus_g_dy2, img_data); //(A+B)^2

    //R = detH - lambda * (TrH)2
    harris = harris_calc(g_dx2, g_dy2, g_dxdy2, g_dx2_plus_g_dy2_square, img_data);

    //R >= max(R) * th;
    //条件を満たす画素をyuvファイルにて色づけて出力
    harris_feature(harris, harris_edge, img_data);
    write_yuv_harris(harris, harris_edge, org, argv[2]);

    free_image(org);

    return 0;
}
