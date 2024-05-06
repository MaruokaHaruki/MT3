#include <Novice.h>
//#include <math.h>
#include <cmath>
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
	// 単位行列を初期化
	Matrix4x4 result = IdentityMatrix();

	// 拡大縮小成分をセット
	result.m4x4[0][0] = scale.x;
	result.m4x4[1][1] = scale.y;
	result.m4x4[2][2] = scale.z;

	return result;
}

//X軸回転行列(yz)
Matrix4x4 MakeRotateXMatrix(float radian) {
	// 単位行列で初期化
	Matrix4x4 result = IdentityMatrix();

	//行列の計算
	result.m4x4[1][1] = std::cos(radian);
	result.m4x4[1][2] = std::sin(radian);
	result.m4x4[2][1] = -std::sin(radian);
	result.m4x4[2][2] = std::cos(radian);

	return result;
}

//Y軸回転行列(zx)
Matrix4x4 MakeRotateYMatrix(float radian) {
	// 単位行列で初期化
	Matrix4x4 result = IdentityMatrix();

	//行列の計算
	result.m4x4[0][0] = std::cos(radian);
	result.m4x4[0][2] = -std::sin(radian);
	result.m4x4[2][0] = std::sin(radian);
	result.m4x4[2][2] = std::cos(radian);

	return result;
}

//Z軸回転行列(xy)
Matrix4x4 MakeRotateZMatrix(float radian) {
	// 単位行列で初期化
	Matrix4x4 result = IdentityMatrix();

	//行列の計算
	result.m4x4[0][0] = std::cos(radian);
	result.m4x4[0][1] = std::sin(radian);
	result.m4x4[1][0] = -std::sin(radian);
	result.m4x4[1][1] = std::cos(radian);

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

//アフィン変換
Matrix4x4 MakeAffineMatrix(const Vector3& scale, const Vector3& rotate, const Vector3& translate) {
	//縮小拡大
	Matrix4x4 scaleMatrix = MakeScaleMatrix(scale);
	//回転
	Matrix4x4 rotateXMatrix = MakeRotateXMatrix(rotate.x);
	Matrix4x4 rotateYMatrix = MakeRotateYMatrix(rotate.y);
	Matrix4x4 rotateZMatrix = MakeRotateZMatrix(rotate.z);
	Matrix4x4 rotateXYZMatrix = MultiplyMatrix(rotateXMatrix, MultiplyMatrix(rotateYMatrix, rotateZMatrix));
	//並行移動
	Matrix4x4 translateMatrix = MakeTranslateMatrix(translate);
	//合成
	Matrix4x4 result = IdentityMatrix();
	result = MultiplyMatrix(scaleMatrix, rotateXYZMatrix);
	result = MultiplyMatrix(result, translateMatrix);

	return result;
}


#pragma endregion

///
///
///

///
///レンダリングパイプライン
///


//正射影行列
Matrix4x4 MakeOrthographicMatrix(float left, float top, float right, float bottom, float nearClip, float farClip) {
	// 単位行列で初期化
	Matrix4x4 result = IdentityMatrix();

	// 正射影平面の範囲から正射影行列を構築する
	result.m4x4[0][0] = 2.0f / ( right - left );
	result.m4x4[1][1] = 2.0f / ( top - bottom );
	result.m4x4[2][2] = -2.0f / ( farClip - nearClip );
	result.m4x4[3][0] = -( right + left ) / ( right - left );
	result.m4x4[3][1] = -( top + bottom ) / ( top - bottom );
	result.m4x4[3][2] = -( farClip + nearClip ) / ( farClip - nearClip );

	return result;

}

//透視投影行列
Matrix4x4 MakePerspectiveFovMatrix(float fovY, float aspectRatio, float nearClip, float farClip) {
	// 単位行列で初期化
	Matrix4x4 result = IdentityMatrix();

	// tanの逆説
	float cot = tanf(fovY / 2.0f);

	// 射影変換行列の要素を計算する
	result.m4x4[0][0] = 1.0f / ( aspectRatio * cot );
	result.m4x4[1][1] = 1.0f / cot;
	result.m4x4[2][2] = farClip / ( farClip - nearClip );
	result.m4x4[2][3] = 1.0f;
	result.m4x4[3][2] = -nearClip * farClip / ( farClip - nearClip );
	result.m4x4[3][3] = 0.0f;

	return result;
}

//ビューポート行列
Matrix4x4 MakeViewportMatrix(float left, float top, float width, float height, float minDepth, float maxDepth) {
	// 単位行列で初期化
	Matrix4x4 result = IdentityMatrix();

	// ビューポート変換行列の要素を計算する
	result.m4x4[0][0] = width / 2.0f;
	result.m4x4[1][1] = -height / 2.0f;
	result.m4x4[2][2] = maxDepth - minDepth;
	result.m4x4[3][0] = left + width / 2.0f;
	result.m4x4[3][1] = top + height / 2.0f;
	result.m4x4[3][2] = minDepth;

	return result;
}

///
///
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

	Matrix4x4 orthographicMatrix = MakeOrthographicMatrix(-160.0f, 160.0f, 200.0f, 300.0f, 0.0f, 1000.0f);

	Matrix4x4 perspectiveFovMatrix = MakePerspectiveFovMatrix(0.63f, 1.33f, 0.1f, 1000.0f);

	Matrix4x4 viewportMatrix = MakeViewportMatrix(100.0f,200.0f,600.0f,300.0f,0.0f,1.0f);


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
		///数値の描画
		///-------------------------------

	/*	Matrix4x4ScreenPrintf(0, 0, worldMatrix, "worldMatrix");*/

		Matrix4x4ScreenPrintf(0, 0, orthographicMatrix, "orthographicMatrix");
		Matrix4x4ScreenPrintf(0, kRowHeight * 5, perspectiveFovMatrix, "perspectiveFovMatrix");
		Matrix4x4ScreenPrintf(0, kRowHeight * 10, viewportMatrix, "viewportMatrix");


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
