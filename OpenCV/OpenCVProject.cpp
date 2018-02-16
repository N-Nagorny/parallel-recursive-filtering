#include <iostream>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

using namespace std;
using namespace cv;

Mat filter(Mat srcImg, int K, int j_k[], double A[], int N, int M_min, int M_max) {
	Mat g = srcImg;
	Mat D33 = Mat::zeros(g.rows, g.cols, CV_64FC3);

	double gr, g1r, g2r, f1r, f2r, f3r, f4r, f1r2, f2r2;
	double gg, g1g, g2g, f1g, f2g, f3g, f4g, f1g2, f2g2;
	double gb, g1b, g2b, f1b, f2b, f3b, f4b, f1b2, f2b2;

	for (int i = 0; i< K; i++)
	{
		double c = 2 * cos(CV_PI*j_k[i] / N);

		for (int y = 0; y <g.rows; y++)
		{
			gr = g1r = g2r = f1r = f2r = f3r = f4r = f1r2 = f2r2 = 0;
			gg = g1g = g2g = f1g = f2g = f3g = f4g = f1g2 = f2g2 = 0;
			gb = g1b = g2b = f1b = f2b = f3b = f4b = f1b2 = f2b2 = 0;

			Vec3b* p = g.ptr<Vec3b>(y);

			for (int x = -M_min; x <g.cols; x++)
			{
				// проверка выхода за пределы слева (всё, что левее, равно нулю)
				if (x - M_max - 1 < 0)
				{
					f2r = f1r = p[x + M_min][0];
					f2g = f1g = p[x + M_min][1];
					f2b = f1b = p[x + M_min][2];
				}
				// проверка выхода за пределы справа (всё, что правее, равно крайнему правому пикселю)
				else if (x + M_min >= g.cols)
				{
					f1r = p[g.cols - 1][0] - p[x - M_max - 1][0];
					f1g = p[g.cols - 1][1] - p[x - M_max - 1][1];
					f1b = p[g.cols - 1][2] - p[x - M_max - 1][2];

					f2r = p[g.cols - 1][0] + p[x - M_max - 1][0];
					f2g = p[g.cols - 1][1] + p[x - M_max - 1][1];
					f2b = p[g.cols - 1][2] + p[x - M_max - 1][2];
				}
				else
				{
					f1r = p[x + M_min][0] - p[x - M_max - 1][0];
					f1g = p[x + M_min][1] - p[x - M_max - 1][1];
					f1b = p[x + M_min][2] - p[x - M_max - 1][2];

					f2r = p[x + M_min][0] + p[x - M_max - 1][0];
					f2g = p[x + M_min][1] + p[x - M_max - 1][1];
					f2b = p[x + M_min][2] + p[x - M_max - 1][2];
				}

				f3r = f1r - f1r2;
				f3g = f1g - f1g2;
				f3b = f1b - f1b2;

				f4r = f2r - f2r2;
				f4g = f2g - f2g2;
				f4b = f2b - f2b2;

				if (j_k[i] == 0)
				{
					gr = g1r + f1r;
					gg = g1g + f1g;
					gb = g1b + f1b;
				}
				//если j_k чётное
				else if ((j_k[i] % 2) == 0)
				{
					gr = c*g1r - g2r + f3r;
					gg = c*g1g - g2g + f3g;
					gb = c*g1b - g2b + f3b;
				}
				//если j_k нечётное
				else
				{
					gr = c*g1r - g2r + f4r;
					gg = c*g1g - g2g + f4g;
					gb = c*g1b - g2b + f4b;
				}

				//складываем результаты обработки звеньев
				if (x >= 0)
				{
					D33.at<Vec3d>(y, x)[0] += gr*A[i];
					D33.at<Vec3d>(y, x)[1] += gg*A[i];
					D33.at<Vec3d>(y, x)[2] += gb*A[i];
				}

				//запоминаем предыдущие расчеты
				g2r = g1r; g1r = gr;
				g2g = g1g; g1g = gg;
				g2b = g1b; g1b = gb;
				f1r2 = f1r; f1g2 = f1g; f1b2 = f1b;
				f2r2 = f2r; f2g2 = f2g; f2b2 = f2b;
			}
		}

	}

	D33.convertTo(g, srcImg.type());
	return g;
}

Mat bpass_filter(Mat srcImg, int N, int M_min, int M_max) {
	const int K = 26; // amount of non-null parallel-recursive filter's cells
	int j_k[K] = { 6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31 };
	double A[K] = {
		-0.0502,
		0.0219,
		0.0557,
		-0.0227,
		-0.0471,
		0.0287,
		0.0444,
		-0.0290,
		-0.0364,
		0.0319,
		0.0320,
		-0.0305,
		-0.0246,
		0.0304,
		0.0197,
		-0.0272,
		-0.0136,
		0.0244,
		0.0094,
		-0.0197,
		-0.0052,
		0.0148,
		0.0024,
		-0.0090,
		-0.0006,
		0.0031 };
	return filter(srcImg, K, j_k, A, N, M_min, M_max);

}

Mat bstop_filter(Mat srcImg, int N, int M_min, int M_max) {
	const int K = 14; // amount of non-null parallel-recursive filter's cells
	int j_k[K];
	for (int i = 0; i < K; ++i) {
		j_k[i] = i + 1;
	}
	double A[K] = {
		0.0301982892696255,
		-0.00305609206311830,
		-0.0640164219031687,
		0.00913565483138058,
		0.0579551788728367,
		-0.0145412099060867,
		-0.0595191875447071,
		0.0202414565600136,
		0.0507651223923958,
		-0.0233147069378204,
		-0.0520974208159894,
		0.0294090190785515,
		0.0336069788856744,
		-0.00951536061584737 };
	return filter(srcImg, K, j_k, A, N, M_min, M_max);
}

int main(int argc, char *argv[])
{
	int filter_number;
	cout << "Write '1' for band-pass filter and '2' for band-stop filter." << endl;
	cin >> filter_number;
	Mat image = imread("image.jpg", 1);
	imshow("figure", image);
	int M_min = 16;
	int M_max = 15;
	int N = M_min + M_max + 1; // amount of samples of Fourier transform of impulse response
	switch (filter_number) {
	case 1:
		imshow("gg", bpass_filter(image, N, M_min, M_max));
		break;
	case 2:
		imshow("gg", bstop_filter(image, N, M_min, M_max));
		break;
	default:
		imshow("gg", bpass_filter(image, N, M_min, M_max));
		break;
	}
	while (waitKey(0) != 27);
	return 0;
}
