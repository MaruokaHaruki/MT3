#include <Novice.h>
//#include <math.h>
#include<Vector3.h>
#include "assert.h"

const char kWindowTitle[] = "LE2B_23_マルオカ_ハルキ";


///-------------------------------
///行列の宣言
///-------------------------------
struct Matrix4x4 {
	float m4x4[4][4];
};

///-------------------------------
///関数の宣言
///-------------------------------

///
///4x4の計算
///

#pragma region 4x4の計算

//1.行列の加法
Matrix4x4 AddMatrix(const Matrix4x4& m1, const Matrix4x4& m2) {
	Matrix4x4 result;
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			result.m4x4[i][j] = m1.m4x4[i][j] + m2.m4x4[i][j];
		}
	}
	return result;
}

//2.行列の減算
Matrix4x4 SubtractMatrix(const Matrix4x4& m1, const Matrix4x4& m2) {
	Matrix4x4 result;
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			result.m4x4[i][j] = m1.m4x4[i][j] - m2.m4x4[i][j];
		}
	}
	return result;
}

//3.行列の積
Matrix4x4 MultiplyMatrix(const Matrix4x4& m1, const Matrix4x4& m2) {
	Matrix4x4 result;

	// 行列の各成分を直接計算して結果を求める
	result.m4x4[0][0] = m1.m4x4[0][0] * m2.m4x4[0][0] + m1.m4x4[0][1] * m2.m4x4[1][0] + m1.m4x4[0][2] * m2.m4x4[2][0] + m1.m4x4[0][3] * m2.m4x4[3][0];
	result.m4x4[0][1] = m1.m4x4[0][0] * m2.m4x4[0][1] + m1.m4x4[0][1] * m2.m4x4[1][1] + m1.m4x4[0][2] * m2.m4x4[2][1] + m1.m4x4[0][3] * m2.m4x4[3][1];
	result.m4x4[0][2] = m1.m4x4[0][0] * m2.m4x4[0][2] + m1.m4x4[0][1] * m2.m4x4[1][2] + m1.m4x4[0][2] * m2.m4x4[2][2] + m1.m4x4[0][3] * m2.m4x4[3][2];
	result.m4x4[0][3] = m1.m4x4[0][0] * m2.m4x4[0][3] + m1.m4x4[0][1] * m2.m4x4[1][3] + m1.m4x4[0][2] * m2.m4x4[2][3] + m1.m4x4[0][3] * m2.m4x4[3][3];

	result.m4x4[1][0] = m1.m4x4[1][0] * m2.m4x4[0][0] + m1.m4x4[1][1] * m2.m4x4[1][0] + m1.m4x4[1][2] * m2.m4x4[2][0] + m1.m4x4[1][3] * m2.m4x4[3][0];
	result.m4x4[1][1] = m1.m4x4[1][0] * m2.m4x4[0][1] + m1.m4x4[1][1] * m2.m4x4[1][1] + m1.m4x4[1][2] * m2.m4x4[2][1] + m1.m4x4[1][3] * m2.m4x4[3][1];
	result.m4x4[1][2] = m1.m4x4[1][0] * m2.m4x4[0][2] + m1.m4x4[1][1] * m2.m4x4[1][2] + m1.m4x4[1][2] * m2.m4x4[2][2] + m1.m4x4[1][3] * m2.m4x4[3][2];
	result.m4x4[1][3] = m1.m4x4[1][0] * m2.m4x4[0][3] + m1.m4x4[1][1] * m2.m4x4[1][3] + m1.m4x4[1][2] * m2.m4x4[2][3] + m1.m4x4[1][3] * m2.m4x4[3][3];

	result.m4x4[2][0] = m1.m4x4[2][0] * m2.m4x4[0][0] + m1.m4x4[2][1] * m2.m4x4[1][0] + m1.m4x4[2][2] * m2.m4x4[2][0] + m1.m4x4[2][3] * m2.m4x4[3][0];
	result.m4x4[2][1] = m1.m4x4[2][0] * m2.m4x4[0][1] + m1.m4x4[2][1] * m2.m4x4[1][1] + m1.m4x4[2][2] * m2.m4x4[2][1] + m1.m4x4[2][3] * m2.m4x4[3][1];
	result.m4x4[2][2] = m1.m4x4[2][0] * m2.m4x4[0][2] + m1.m4x4[2][1] * m2.m4x4[1][2] + m1.m4x4[2][2] * m2.m4x4[2][2] + m1.m4x4[2][3] * m2.m4x4[3][2];
	result.m4x4[2][3] = m1.m4x4[2][0] * m2.m4x4[0][3] + m1.m4x4[2][1] * m2.m4x4[1][3] + m1.m4x4[2][2] * m2.m4x4[2][3] + m1.m4x4[2][3] * m2.m4x4[3][3];

	result.m4x4[3][0] = m1.m4x4[3][0] * m2.m4x4[0][0] + m1.m4x4[3][1] * m2.m4x4[1][0] + m1.m4x4[3][2] * m2.m4x4[2][0] + m1.m4x4[3][3] * m2.m4x4[3][0];
	result.m4x4[3][1] = m1.m4x4[3][0] * m2.m4x4[0][1] + m1.m4x4[3][1] * m2.m4x4[1][1] + m1.m4x4[3][2] * m2.m4x4[2][1] + m1.m4x4[3][3] * m2.m4x4[3][1];
	result.m4x4[3][2] = m1.m4x4[3][0] * m2.m4x4[0][2] + m1.m4x4[3][1] * m2.m4x4[1][2] + m1.m4x4[3][2] * m2.m4x4[2][2] + m1.m4x4[3][3] * m2.m4x4[3][2];
	result.m4x4[3][3] = m1.m4x4[3][0] * m2.m4x4[0][3] + m1.m4x4[3][1] * m2.m4x4[1][3] + m1.m4x4[3][2] * m2.m4x4[2][3] + m1.m4x4[3][3] * m2.m4x4[3][3];

	return result;
}

// 逆行列を計算する関数
Matrix4x4 InverseMatrix(const Matrix4x4& matrix) {
	// 行列式を計算
	float det =
		matrix.m4x4[0][0] * ( matrix.m4x4[1][1] * matrix.m4x4[2][2] * matrix.m4x4[3][3] +
			matrix.m4x4[1][2] * matrix.m4x4[2][3] * matrix.m4x4[3][1] +
			matrix.m4x4[1][3] * matrix.m4x4[2][1] * matrix.m4x4[3][2] -
			matrix.m4x4[1][3] * matrix.m4x4[2][2] * matrix.m4x4[3][1] -
			matrix.m4x4[1][1] * matrix.m4x4[2][3] * matrix.m4x4[3][2] -
			matrix.m4x4[1][2] * matrix.m4x4[2][1] * matrix.m4x4[3][3] ) -
		matrix.m4x4[0][1] * ( matrix.m4x4[1][0] * matrix.m4x4[2][2] * matrix.m4x4[3][3] +
			matrix.m4x4[1][2] * matrix.m4x4[2][3] * matrix.m4x4[3][0] +
			matrix.m4x4[1][3] * matrix.m4x4[2][0] * matrix.m4x4[3][2] -
			matrix.m4x4[1][3] * matrix.m4x4[2][2] * matrix.m4x4[3][0] -
			matrix.m4x4[1][0] * matrix.m4x4[2][3] * matrix.m4x4[3][2] -
			matrix.m4x4[1][2] * matrix.m4x4[2][0] * matrix.m4x4[3][3] ) +
		matrix.m4x4[0][2] * ( matrix.m4x4[1][0] * matrix.m4x4[2][1] * matrix.m4x4[3][3] +
			matrix.m4x4[1][1] * matrix.m4x4[2][3] * matrix.m4x4[3][0] +
			matrix.m4x4[1][3] * matrix.m4x4[2][0] * matrix.m4x4[3][1] -
			matrix.m4x4[1][3] * matrix.m4x4[2][1] * matrix.m4x4[3][0] -
			matrix.m4x4[1][0] * matrix.m4x4[2][3] * matrix.m4x4[3][1] -
			matrix.m4x4[1][1] * matrix.m4x4[2][0] * matrix.m4x4[3][3] ) -
		matrix.m4x4[0][3] * ( matrix.m4x4[1][0] * matrix.m4x4[2][1] * matrix.m4x4[3][2] +
			matrix.m4x4[1][1] * matrix.m4x4[2][2] * matrix.m4x4[3][0] +
			matrix.m4x4[1][2] * matrix.m4x4[2][0] * matrix.m4x4[3][1] -
			matrix.m4x4[1][2] * matrix.m4x4[2][1] * matrix.m4x4[3][0] -
			matrix.m4x4[1][0] * matrix.m4x4[2][2] * matrix.m4x4[3][1] -
			matrix.m4x4[1][1] * matrix.m4x4[2][0] * matrix.m4x4[3][2] );

	Matrix4x4 result;
	// 各要素について余因子を計算して逆行列を出力
	result.m4x4[0][0] = ( matrix.m4x4[1][1] * matrix.m4x4[2][2] * matrix.m4x4[3][3] +
		matrix.m4x4[1][2] * matrix.m4x4[2][3] * matrix.m4x4[3][1] +
		matrix.m4x4[1][3] * matrix.m4x4[2][1] * matrix.m4x4[3][2] -
		matrix.m4x4[1][3] * matrix.m4x4[2][2] * matrix.m4x4[3][1] -
		matrix.m4x4[1][1] * matrix.m4x4[2][3] * matrix.m4x4[3][2] -
		matrix.m4x4[1][2] * matrix.m4x4[2][1] * matrix.m4x4[3][3] ) / det;

	result.m4x4[0][1] = ( matrix.m4x4[0][1] * matrix.m4x4[2][3] * matrix.m4x4[3][2] +
		matrix.m4x4[0][2] * matrix.m4x4[2][1] * matrix.m4x4[3][3] +
		matrix.m4x4[0][3] * matrix.m4x4[2][2] * matrix.m4x4[3][1] -
		matrix.m4x4[0][3] * matrix.m4x4[2][1] * matrix.m4x4[3][2] -
		matrix.m4x4[0][1] * matrix.m4x4[2][2] * matrix.m4x4[3][3] -
		matrix.m4x4[0][2] * matrix.m4x4[2][3] * matrix.m4x4[3][1] ) / det;

	result.m4x4[0][2] = ( matrix.m4x4[0][1] * matrix.m4x4[1][2] * matrix.m4x4[3][3] +
		matrix.m4x4[0][2] * matrix.m4x4[1][3] * matrix.m4x4[3][1] +
		matrix.m4x4[0][3] * matrix.m4x4[1][1] * matrix.m4x4[3][2] -
		matrix.m4x4[0][3] * matrix.m4x4[1][2] * matrix.m4x4[3][1] -
		matrix.m4x4[0][1] * matrix.m4x4[1][3] * matrix.m4x4[3][2] -
		matrix.m4x4[0][2] * matrix.m4x4[1][1] * matrix.m4x4[3][3] ) / det;

	result.m4x4[0][3] = ( matrix.m4x4[0][1] * matrix.m4x4[1][3] * matrix.m4x4[2][2] +
		matrix.m4x4[0][2] * matrix.m4x4[1][1] * matrix.m4x4[2][3] +
		matrix.m4x4[0][3] * matrix.m4x4[1][2] * matrix.m4x4[2][1] -
		matrix.m4x4[0][3] * matrix.m4x4[1][1] * matrix.m4x4[2][2] -
		matrix.m4x4[0][1] * matrix.m4x4[1][2] * matrix.m4x4[2][3] -
		matrix.m4x4[0][2] * matrix.m4x4[1][3] * matrix.m4x4[2][1] ) / det;

	result.m4x4[1][0] = ( matrix.m4x4[1][0] * matrix.m4x4[2][2] * matrix.m4x4[3][3] +
		matrix.m4x4[1][2] * matrix.m4x4[2][3] * matrix.m4x4[3][0] +
		matrix.m4x4[1][3] * matrix.m4x4[2][0] * matrix.m4x4[3][2] -
		matrix.m4x4[1][3] * matrix.m4x4[2][2] * matrix.m4x4[3][0] -
		matrix.m4x4[1][0] * matrix.m4x4[2][3] * matrix.m4x4[3][2] -
		matrix.m4x4[1][2] * matrix.m4x4[2][0] * matrix.m4x4[3][3] ) / det;

	result.m4x4[1][1] = ( matrix.m4x4[0][0] * matrix.m4x4[2][3] * matrix.m4x4[3][2] +
		matrix.m4x4[0][2] * matrix.m4x4[2][0] * matrix.m4x4[3][3] +
		matrix.m4x4[0][3] * matrix.m4x4[2][2] * matrix.m4x4[3][0] -
		matrix.m4x4[0][3] * matrix.m4x4[2][0] * matrix.m4x4[3][2] -
		matrix.m4x4[0][0] * matrix.m4x4[2][2] * matrix.m4x4[3][3] -
		matrix.m4x4[0][2] * matrix.m4x4[2][3] * matrix.m4x4[3][0] ) / det;

	result.m4x4[1][2] = ( matrix.m4x4[0][0] * matrix.m4x4[1][2] * matrix.m4x4[3][3] +
		matrix.m4x4[0][2] * matrix.m4x4[1][3] * matrix.m4x4[3][0] +
		matrix.m4x4[0][3] * matrix.m4x4[1][0] * matrix.m4x4[3][2] -
		matrix.m4x4[0][3] * matrix.m4x4[1][2] * matrix.m4x4[3][0] -
		matrix.m4x4[0][0] * matrix.m4x4[1][3] * matrix.m4x4[3][2] -
		matrix.m4x4[0][2] * matrix.m4x4[1][0] * matrix.m4x4[3][3] ) / det;

	result.m4x4[1][3] = ( matrix.m4x4[0][0] * matrix.m4x4[1][3] * matrix.m4x4[2][2] +
		matrix.m4x4[0][2] * matrix.m4x4[1][0] * matrix.m4x4[2][3] +
		matrix.m4x4[0][3] * matrix.m4x4[1][2] * matrix.m4x4[2][0] -
		matrix.m4x4[0][3] * matrix.m4x4[1][0] * matrix.m4x4[2][2] -
		matrix.m4x4[0][0] * matrix.m4x4[1][2] * matrix.m4x4[2][3] -
		matrix.m4x4[0][2] * matrix.m4x4[1][3] * matrix.m4x4[2][0] ) / det;

	result.m4x4[2][0] = ( matrix.m4x4[1][0] * matrix.m4x4[2][1] * matrix.m4x4[3][3] +
		matrix.m4x4[1][1] * matrix.m4x4[2][3] * matrix.m4x4[3][0] +
		matrix.m4x4[1][3] * matrix.m4x4[2][0] * matrix.m4x4[3][1] -
		matrix.m4x4[1][3] * matrix.m4x4[2][1] * matrix.m4x4[3][0] -
		matrix.m4x4[1][0] * matrix.m4x4[2][3] * matrix.m4x4[3][1] -
		matrix.m4x4[1][1] * matrix.m4x4[2][0] * matrix.m4x4[3][3] ) / det;

	result.m4x4[2][1] = ( matrix.m4x4[0][0] * matrix.m4x4[2][3] * matrix.m4x4[3][1] +
		matrix.m4x4[0][1] * matrix.m4x4[2][0] * matrix.m4x4[3][3] +
		matrix.m4x4[0][3] * matrix.m4x4[2][1] * matrix.m4x4[3][0] -
		matrix.m4x4[0][3] * matrix.m4x4[2][0] * matrix.m4x4[3][1] -
		matrix.m4x4[0][0] * matrix.m4x4[2][1] * matrix.m4x4[3][3] -
		matrix.m4x4[0][1] * matrix.m4x4[2][3] * matrix.m4x4[3][0] ) / det;

	result.m4x4[2][2] = ( matrix.m4x4[0][0] * matrix.m4x4[1][1] * matrix.m4x4[3][3] +
		matrix.m4x4[0][1] * matrix.m4x4[1][3] * matrix.m4x4[3][0] +
		matrix.m4x4[0][3] * matrix.m4x4[1][0] * matrix.m4x4[3][1] -
		matrix.m4x4[0][3] * matrix.m4x4[1][1] * matrix.m4x4[3][0] -
		matrix.m4x4[0][0] * matrix.m4x4[1][3] * matrix.m4x4[3][1] -
		matrix.m4x4[0][1] * matrix.m4x4[1][0] * matrix.m4x4[3][3] ) / det;

	result.m4x4[2][3] = ( matrix.m4x4[0][0] * matrix.m4x4[1][3] * matrix.m4x4[2][1] +
		matrix.m4x4[0][1] * matrix.m4x4[1][0] * matrix.m4x4[2][3] +
		matrix.m4x4[0][3] * matrix.m4x4[1][1] * matrix.m4x4[2][0] -
		matrix.m4x4[0][3] * matrix.m4x4[1][0] * matrix.m4x4[2][1] -
		matrix.m4x4[0][0] * matrix.m4x4[1][1] * matrix.m4x4[2][3] -
		matrix.m4x4[0][1] * matrix.m4x4[1][3] * matrix.m4x4[2][0] ) / det;

	result.m4x4[3][0] = ( matrix.m4x4[1][0] * matrix.m4x4[2][2] * matrix.m4x4[3][1] +
		matrix.m4x4[1][1] * matrix.m4x4[2][0] * matrix.m4x4[3][2] +
		matrix.m4x4[1][2] * matrix.m4x4[2][1] * matrix.m4x4[3][0] -
		matrix.m4x4[1][2] * matrix.m4x4[2][0] * matrix.m4x4[3][1] -
		matrix.m4x4[1][0] * matrix.m4x4[2][1] * matrix.m4x4[3][2] -
		matrix.m4x4[1][1] * matrix.m4x4[2][2] * matrix.m4x4[3][0] ) / det;

	result.m4x4[3][1] = ( matrix.m4x4[0][0] * matrix.m4x4[2][1] * matrix.m4x4[3][2] +
		matrix.m4x4[0][1] * matrix.m4x4[2][2] * matrix.m4x4[3][0] +
		matrix.m4x4[0][2] * matrix.m4x4[2][0] * matrix.m4x4[3][1] -
		matrix.m4x4[0][2] * matrix.m4x4[2][1] * matrix.m4x4[3][0] -
		matrix.m4x4[0][0] * matrix.m4x4[2][2] * matrix.m4x4[3][1] -
		matrix.m4x4[0][1] * matrix.m4x4[2][0] * matrix.m4x4[3][2] ) / det;

	result.m4x4[3][2] = ( matrix.m4x4[0][0] * matrix.m4x4[1][2] * matrix.m4x4[3][1] +
		matrix.m4x4[0][1] * matrix.m4x4[1][0] * matrix.m4x4[3][2] +
		matrix.m4x4[0][2] * matrix.m4x4[1][1] * matrix.m4x4[3][0] -
		matrix.m4x4[0][2] * matrix.m4x4[1][0] * matrix.m4x4[3][1] -
		matrix.m4x4[0][0] * matrix.m4x4[1][1] * matrix.m4x4[3][2] -
		matrix.m4x4[0][1] * matrix.m4x4[1][2] * matrix.m4x4[3][0] ) / det;

	result.m4x4[3][3] = ( matrix.m4x4[0][0] * matrix.m4x4[1][1] * matrix.m4x4[2][2] +
		matrix.m4x4[0][1] * matrix.m4x4[1][2] * matrix.m4x4[2][0] +
		matrix.m4x4[0][2] * matrix.m4x4[1][0] * matrix.m4x4[2][1] -
		matrix.m4x4[0][2] * matrix.m4x4[1][1] * matrix.m4x4[2][0] -
		matrix.m4x4[0][0] * matrix.m4x4[1][2] * matrix.m4x4[2][1] -
		matrix.m4x4[0][1] * matrix.m4x4[1][0] * matrix.m4x4[2][2] ) / det;

	return result;
}

//5.転置行列
Matrix4x4 TransposeMatrix(const Matrix4x4& matrix) {
	Matrix4x4 result;
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			result.m4x4[i][j] = matrix.m4x4[j][i];
		}
	}
	return result;
}

//6.単位行列の作成
Matrix4x4 IdentityMatrix() {
	Matrix4x4 result;
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			if (i == j) {
				result.m4x4[i][j] = 1.0f;
			} else {
				result.m4x4[i][j] = 0.0f;
			}
		}
	}
	return result;
}

#pragma endregion

///
///end
///


///
///3次元アフィン行列
///

#pragma region 3次元アフィン行列


// 3次元ベクトル表示
static void Vector3ScreenPrintf(int x, int y, Vector3 vector3, const char* label) {
	Novice::ScreenPrintf(x, y, "%4.2f", vector3.x);
	Novice::ScreenPrintf(x + 50, y, "%4.2f", vector3.y);
	Novice::ScreenPrintf(x + 100, y, "%4.2f", vector3.z);
	Novice::ScreenPrintf(x + 150, y, "%s", label);
};

// 4x4行列の数値表示
static const int kRowHeight = 20;
static const int kColumnWight = 60;
void Matrix4x4ScreenPrintf(int x, int y, const Matrix4x4& matrix, const char* label) {
	for (int row = 0; row < 4; ++row) {
		for (int column = 0; column < 4; ++column) {
			Novice::ScreenPrintf(x + column * kColumnWight, y + row * kRowHeight, "%6.02f", matrix.m4x4[row][column]);
		}
	}

	Novice::ScreenPrintf(x + 250, y, "%s", label);
}

// 単位行列の作成
Matrix4x4 IdentityMatrix() {
	Matrix4x4 result;
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			if (i == j) {
				result.m4x4[i][j] = 1.0f;
			} else {
				result.m4x4[i][j] = 0.0f;
			}
		}
	}
	return result;
}

// 平行移動行列を作成する関数
Matrix4x4 MakeTranslateMatrix(const Vector3& translate) {
	Matrix4x4 result = IdentityMatrix(); // 単位行列を初期化

	// 平行移動成分をセット
	result.m4x4[3][0] = translate.x;
	result.m4x4[3][1] = translate.y;
	result.m4x4[3][2] = translate.z;

	return result;
}

// 拡大縮小行列を作成する関数
Matrix4x4 MakeScaleMatrix(const Vector3& scale) {
	Matrix4x4 result = IdentityMatrix(); // 単位行列を初期化

	// 拡大縮小成分をセット
	result.m4x4[0][0] = scale.x;
	result.m4x4[1][1] = scale.y;
	result.m4x4[2][2] = scale.z;

	return result;
}

// 座標変換を行う関数
Vector3 Transform(const Vector3& vector, const Matrix4x4& matrix) {
	Vector3 result;

	// 行列とベクトルの乗算
	result.x = vector.x * matrix.m4x4[0][0] + vector.y * matrix.m4x4[1][0] + vector.z * matrix.m4x4[2][0] + matrix.m4x4[3][0];
	result.y = vector.x * matrix.m4x4[0][1] + vector.y * matrix.m4x4[1][1] + vector.z * matrix.m4x4[2][1] + matrix.m4x4[3][1];
	result.z = vector.x * matrix.m4x4[0][2] + vector.y * matrix.m4x4[1][2] + vector.z * matrix.m4x4[2][2] + matrix.m4x4[3][2];
	float w = vector.x * matrix.m4x4[0][3] + vector.y * matrix.m4x4[1][3] + vector.z * matrix.m4x4[2][3] + matrix.m4x4[3][3];
	//この処理を忘れない！
	assert(w != 0.0f);
	result.x /= w;
	result.y /= w;
	result.z /= w;
	return result;
}

#pragma endregion

///
///end
///


// Windowsアプリでのエントリーポイント(main関数)
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {

	// ライブラリの初期化
	Novice::Initialize(kWindowTitle, 1280, 720);

	// キー入力結果を受け取る箱
	char keys[256] = { 0 };
	char preKeys[256] = { 0 };

	///-------------------------------
	///変数の宣言
	///-------------------------------
	Vector3 translate{ 4.1f,2.6f,0.8f };
	Vector3 scale{ 1.5f,5.2f,7.3f };
	Matrix4x4 translateMatrix = MakeTranslateMatrix(translate);
	Matrix4x4 scaleMatrix = MakeScaleMatrix(scale);
	Vector3 point{ 2.3f,3.8f,1.4f };
	Matrix4x4 transfromMatrix = {
		1.0f,2.0f,3.0f,4.0f,
		3.0f,1.0f,1.0f,2.0f,
		1.0f,4.0f,2.0f,3.0f,
		2.0f,2.0f,1.0f,3.0f
	};

	Vector3 transformed = Transform(point, transfromMatrix);

	// ウィンドウの×ボタンが押されるまでループ
	while (Novice::ProcessMessage() == 0) {
		// フレームの開始
		Novice::BeginFrame();

		// キー入力を受け取る
		memcpy(preKeys, keys, 256);
		Novice::GetHitKeyStateAll(keys);

		///
		/// ↓更新処理ここから
		///

		///
		/// ↑更新処理ここまで
		///

		///
		/// ↓描画処理ここから
		///


		///-------------------------------
		///数値の秒yが
		///-------------------------------
		Vector3ScreenPrintf(0, 20, transformed, "transformed");
		Matrix4x4ScreenPrintf(0, kRowHeight * 5, translateMatrix, "translateMatrix");
		Matrix4x4ScreenPrintf(0, kRowHeight * 5 * 2, scaleMatrix, "scaleMatrix");

		///
		/// ↑描画処理ここまで
		///

		// フレームの終了
		Novice::EndFrame();

		// ESCキーが押されたらループを抜ける
		if (preKeys[DIK_ESCAPE] == 0 && keys[DIK_ESCAPE] != 0) {
			break;
		}
	}

	// ライブラリの終了
	Novice::Finalize();
	return 0;
}
