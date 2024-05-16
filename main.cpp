#include <Novice.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <cmath>
#include<Vector3.h>
#include "assert.h"
#include "Matrix4x4.h"

#include<imgui.h>

const char kWindowTitle[] = "LE2B_23_マルオカ_ハルキ";


///-------------------------------
///行列の宣言
///-------------------------------
// 球
struct Sphere {
	Vector3 center;
	float radius;
};

///-------------------------------
///関数の宣言
///-------------------------------
//ベクター3描画
static const int kColumnWidth = 60;
void Vector3ScreenPrintf(uint32_t x, uint32_t y, Vector3 v, const char* label) {
	Novice::ScreenPrintf(x, y, "%.02f", v.x);
	Novice::ScreenPrintf(x + kColumnWidth, y, "%.02f", v.y);
	Novice::ScreenPrintf(x + kColumnWidth * 2, y, "%.02f", v.z);
	Novice::ScreenPrintf(x + kColumnWidth * 3, y, "%s", label);
}

//内積
float Dot(const Vector3& v1, const Vector3& v2) {
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

// ベクトル減算
Vector3 SubtractVector3(const Vector3& v1, const Vector3& v2) {
	Vector3 result;
	result.x = v1.x - v2.x;
	result.y = v1.y - v2.y;
	result.z = v1.z - v2.z;
	return result;
}

// スカラー乗算
Vector3 MultiplyVector3(const Vector3& v, float scalar) {
	return Vector3{ v.x * scalar, v.y * scalar, v.z * scalar };
}

// ベクトル加算
Vector3 AddVector3(const Vector3& v1, const Vector3& v2) {
	return Vector3{ v1.x + v2.x, v1.y + v2.y, v1.z + v2.z };
}

///
///4x4の計算
///

#pragma region 4x4の計算

//1.行列の加法
Matrix4x4 AddMatrix(const Matrix4x4& m1, const Matrix4x4& m2) {
	Matrix4x4 result;
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			result.m[i][j] = m1.m[i][j] + m2.m[i][j];
		}
	}
	return result;
}

//2.行列の減算
Matrix4x4 SubtractMatrix(const Matrix4x4& m1, const Matrix4x4& m2) {
	Matrix4x4 result;
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			result.m[i][j] = m1.m[i][j] - m2.m[i][j];
		}
	}
	return result;
}

//3.行列の積
Matrix4x4 MultiplyMatrix(const Matrix4x4& m1, const Matrix4x4& m2) {
	Matrix4x4 result;

	// 行列の各成分を直接計算して結果を求める
	result.m[0][0] = m1.m[0][0] * m2.m[0][0] + m1.m[0][1] * m2.m[1][0] + m1.m[0][2] * m2.m[2][0] + m1.m[0][3] * m2.m[3][0];
	result.m[0][1] = m1.m[0][0] * m2.m[0][1] + m1.m[0][1] * m2.m[1][1] + m1.m[0][2] * m2.m[2][1] + m1.m[0][3] * m2.m[3][1];
	result.m[0][2] = m1.m[0][0] * m2.m[0][2] + m1.m[0][1] * m2.m[1][2] + m1.m[0][2] * m2.m[2][2] + m1.m[0][3] * m2.m[3][2];
	result.m[0][3] = m1.m[0][0] * m2.m[0][3] + m1.m[0][1] * m2.m[1][3] + m1.m[0][2] * m2.m[2][3] + m1.m[0][3] * m2.m[3][3];

	result.m[1][0] = m1.m[1][0] * m2.m[0][0] + m1.m[1][1] * m2.m[1][0] + m1.m[1][2] * m2.m[2][0] + m1.m[1][3] * m2.m[3][0];
	result.m[1][1] = m1.m[1][0] * m2.m[0][1] + m1.m[1][1] * m2.m[1][1] + m1.m[1][2] * m2.m[2][1] + m1.m[1][3] * m2.m[3][1];
	result.m[1][2] = m1.m[1][0] * m2.m[0][2] + m1.m[1][1] * m2.m[1][2] + m1.m[1][2] * m2.m[2][2] + m1.m[1][3] * m2.m[3][2];
	result.m[1][3] = m1.m[1][0] * m2.m[0][3] + m1.m[1][1] * m2.m[1][3] + m1.m[1][2] * m2.m[2][3] + m1.m[1][3] * m2.m[3][3];

	result.m[2][0] = m1.m[2][0] * m2.m[0][0] + m1.m[2][1] * m2.m[1][0] + m1.m[2][2] * m2.m[2][0] + m1.m[2][3] * m2.m[3][0];
	result.m[2][1] = m1.m[2][0] * m2.m[0][1] + m1.m[2][1] * m2.m[1][1] + m1.m[2][2] * m2.m[2][1] + m1.m[2][3] * m2.m[3][1];
	result.m[2][2] = m1.m[2][0] * m2.m[0][2] + m1.m[2][1] * m2.m[1][2] + m1.m[2][2] * m2.m[2][2] + m1.m[2][3] * m2.m[3][2];
	result.m[2][3] = m1.m[2][0] * m2.m[0][3] + m1.m[2][1] * m2.m[1][3] + m1.m[2][2] * m2.m[2][3] + m1.m[2][3] * m2.m[3][3];

	result.m[3][0] = m1.m[3][0] * m2.m[0][0] + m1.m[3][1] * m2.m[1][0] + m1.m[3][2] * m2.m[2][0] + m1.m[3][3] * m2.m[3][0];
	result.m[3][1] = m1.m[3][0] * m2.m[0][1] + m1.m[3][1] * m2.m[1][1] + m1.m[3][2] * m2.m[2][1] + m1.m[3][3] * m2.m[3][1];
	result.m[3][2] = m1.m[3][0] * m2.m[0][2] + m1.m[3][1] * m2.m[1][2] + m1.m[3][2] * m2.m[2][2] + m1.m[3][3] * m2.m[3][2];
	result.m[3][3] = m1.m[3][0] * m2.m[0][3] + m1.m[3][1] * m2.m[1][3] + m1.m[3][2] * m2.m[2][3] + m1.m[3][3] * m2.m[3][3];

	return result;
}

// 逆行列を計算する関数
Matrix4x4 InverseMatrix(const Matrix4x4& matrix) {
	// 行列式を計算
	float det =
		matrix.m[0][0] * ( matrix.m[1][1] * matrix.m[2][2] * matrix.m[3][3] +
			matrix.m[1][2] * matrix.m[2][3] * matrix.m[3][1] +
			matrix.m[1][3] * matrix.m[2][1] * matrix.m[3][2] -
			matrix.m[1][3] * matrix.m[2][2] * matrix.m[3][1] -
			matrix.m[1][1] * matrix.m[2][3] * matrix.m[3][2] -
			matrix.m[1][2] * matrix.m[2][1] * matrix.m[3][3] ) -
		matrix.m[0][1] * ( matrix.m[1][0] * matrix.m[2][2] * matrix.m[3][3] +
			matrix.m[1][2] * matrix.m[2][3] * matrix.m[3][0] +
			matrix.m[1][3] * matrix.m[2][0] * matrix.m[3][2] -
			matrix.m[1][3] * matrix.m[2][2] * matrix.m[3][0] -
			matrix.m[1][0] * matrix.m[2][3] * matrix.m[3][2] -
			matrix.m[1][2] * matrix.m[2][0] * matrix.m[3][3] ) +
		matrix.m[0][2] * ( matrix.m[1][0] * matrix.m[2][1] * matrix.m[3][3] +
			matrix.m[1][1] * matrix.m[2][3] * matrix.m[3][0] +
			matrix.m[1][3] * matrix.m[2][0] * matrix.m[3][1] -
			matrix.m[1][3] * matrix.m[2][1] * matrix.m[3][0] -
			matrix.m[1][0] * matrix.m[2][3] * matrix.m[3][1] -
			matrix.m[1][1] * matrix.m[2][0] * matrix.m[3][3] ) -
		matrix.m[0][3] * ( matrix.m[1][0] * matrix.m[2][1] * matrix.m[3][2] +
			matrix.m[1][1] * matrix.m[2][2] * matrix.m[3][0] +
			matrix.m[1][2] * matrix.m[2][0] * matrix.m[3][1] -
			matrix.m[1][2] * matrix.m[2][1] * matrix.m[3][0] -
			matrix.m[1][0] * matrix.m[2][2] * matrix.m[3][1] -
			matrix.m[1][1] * matrix.m[2][0] * matrix.m[3][2] );

	Matrix4x4 result;
	// 各要素について余因子を計算して逆行列を出力
	result.m[0][0] = ( matrix.m[1][1] * matrix.m[2][2] * matrix.m[3][3] +
		matrix.m[1][2] * matrix.m[2][3] * matrix.m[3][1] +
		matrix.m[1][3] * matrix.m[2][1] * matrix.m[3][2] -
		matrix.m[1][3] * matrix.m[2][2] * matrix.m[3][1] -
		matrix.m[1][1] * matrix.m[2][3] * matrix.m[3][2] -
		matrix.m[1][2] * matrix.m[2][1] * matrix.m[3][3] ) / det;

	result.m[0][1] = ( matrix.m[0][1] * matrix.m[2][3] * matrix.m[3][2] +
		matrix.m[0][2] * matrix.m[2][1] * matrix.m[3][3] +
		matrix.m[0][3] * matrix.m[2][2] * matrix.m[3][1] -
		matrix.m[0][3] * matrix.m[2][1] * matrix.m[3][2] -
		matrix.m[0][1] * matrix.m[2][2] * matrix.m[3][3] -
		matrix.m[0][2] * matrix.m[2][3] * matrix.m[3][1] ) / det;

	result.m[0][2] = ( matrix.m[0][1] * matrix.m[1][2] * matrix.m[3][3] +
		matrix.m[0][2] * matrix.m[1][3] * matrix.m[3][1] +
		matrix.m[0][3] * matrix.m[1][1] * matrix.m[3][2] -
		matrix.m[0][3] * matrix.m[1][2] * matrix.m[3][1] -
		matrix.m[0][1] * matrix.m[1][3] * matrix.m[3][2] -
		matrix.m[0][2] * matrix.m[1][1] * matrix.m[3][3] ) / det;

	result.m[0][3] = ( matrix.m[0][1] * matrix.m[1][3] * matrix.m[2][2] +
		matrix.m[0][2] * matrix.m[1][1] * matrix.m[2][3] +
		matrix.m[0][3] * matrix.m[1][2] * matrix.m[2][1] -
		matrix.m[0][3] * matrix.m[1][1] * matrix.m[2][2] -
		matrix.m[0][1] * matrix.m[1][2] * matrix.m[2][3] -
		matrix.m[0][2] * matrix.m[1][3] * matrix.m[2][1] ) / det;

	result.m[1][0] = ( matrix.m[1][0] * matrix.m[2][2] * matrix.m[3][3] +
		matrix.m[1][2] * matrix.m[2][3] * matrix.m[3][0] +
		matrix.m[1][3] * matrix.m[2][0] * matrix.m[3][2] -
		matrix.m[1][3] * matrix.m[2][2] * matrix.m[3][0] -
		matrix.m[1][0] * matrix.m[2][3] * matrix.m[3][2] -
		matrix.m[1][2] * matrix.m[2][0] * matrix.m[3][3] ) / det;

	result.m[1][1] = ( matrix.m[0][0] * matrix.m[2][3] * matrix.m[3][2] +
		matrix.m[0][2] * matrix.m[2][0] * matrix.m[3][3] +
		matrix.m[0][3] * matrix.m[2][2] * matrix.m[3][0] -
		matrix.m[0][3] * matrix.m[2][0] * matrix.m[3][2] -
		matrix.m[0][0] * matrix.m[2][2] * matrix.m[3][3] -
		matrix.m[0][2] * matrix.m[2][3] * matrix.m[3][0] ) / det;

	result.m[1][2] = ( matrix.m[0][0] * matrix.m[1][2] * matrix.m[3][3] +
		matrix.m[0][2] * matrix.m[1][3] * matrix.m[3][0] +
		matrix.m[0][3] * matrix.m[1][0] * matrix.m[3][2] -
		matrix.m[0][3] * matrix.m[1][2] * matrix.m[3][0] -
		matrix.m[0][0] * matrix.m[1][3] * matrix.m[3][2] -
		matrix.m[0][2] * matrix.m[1][0] * matrix.m[3][3] ) / det;

	result.m[1][3] = ( matrix.m[0][0] * matrix.m[1][3] * matrix.m[2][2] +
		matrix.m[0][2] * matrix.m[1][0] * matrix.m[2][3] +
		matrix.m[0][3] * matrix.m[1][2] * matrix.m[2][0] -
		matrix.m[0][3] * matrix.m[1][0] * matrix.m[2][2] -
		matrix.m[0][0] * matrix.m[1][2] * matrix.m[2][3] -
		matrix.m[0][2] * matrix.m[1][3] * matrix.m[2][0] ) / det;

	result.m[2][0] = ( matrix.m[1][0] * matrix.m[2][1] * matrix.m[3][3] +
		matrix.m[1][1] * matrix.m[2][3] * matrix.m[3][0] +
		matrix.m[1][3] * matrix.m[2][0] * matrix.m[3][1] -
		matrix.m[1][3] * matrix.m[2][1] * matrix.m[3][0] -
		matrix.m[1][0] * matrix.m[2][3] * matrix.m[3][1] -
		matrix.m[1][1] * matrix.m[2][0] * matrix.m[3][3] ) / det;

	result.m[2][1] = ( matrix.m[0][0] * matrix.m[2][3] * matrix.m[3][1] +
		matrix.m[0][1] * matrix.m[2][0] * matrix.m[3][3] +
		matrix.m[0][3] * matrix.m[2][1] * matrix.m[3][0] -
		matrix.m[0][3] * matrix.m[2][0] * matrix.m[3][1] -
		matrix.m[0][0] * matrix.m[2][1] * matrix.m[3][3] -
		matrix.m[0][1] * matrix.m[2][3] * matrix.m[3][0] ) / det;

	result.m[2][2] = ( matrix.m[0][0] * matrix.m[1][1] * matrix.m[3][3] +
		matrix.m[0][1] * matrix.m[1][3] * matrix.m[3][0] +
		matrix.m[0][3] * matrix.m[1][0] * matrix.m[3][1] -
		matrix.m[0][3] * matrix.m[1][1] * matrix.m[3][0] -
		matrix.m[0][0] * matrix.m[1][3] * matrix.m[3][1] -
		matrix.m[0][1] * matrix.m[1][0] * matrix.m[3][3] ) / det;

	result.m[2][3] = ( matrix.m[0][0] * matrix.m[1][3] * matrix.m[2][1] +
		matrix.m[0][1] * matrix.m[1][0] * matrix.m[2][3] +
		matrix.m[0][3] * matrix.m[1][1] * matrix.m[2][0] -
		matrix.m[0][3] * matrix.m[1][0] * matrix.m[2][1] -
		matrix.m[0][0] * matrix.m[1][1] * matrix.m[2][3] -
		matrix.m[0][1] * matrix.m[1][3] * matrix.m[2][0] ) / det;

	result.m[3][0] = ( matrix.m[1][0] * matrix.m[2][2] * matrix.m[3][1] +
		matrix.m[1][1] * matrix.m[2][0] * matrix.m[3][2] +
		matrix.m[1][2] * matrix.m[2][1] * matrix.m[3][0] -
		matrix.m[1][2] * matrix.m[2][0] * matrix.m[3][1] -
		matrix.m[1][0] * matrix.m[2][1] * matrix.m[3][2] -
		matrix.m[1][1] * matrix.m[2][2] * matrix.m[3][0] ) / det;

	result.m[3][1] = ( matrix.m[0][0] * matrix.m[2][1] * matrix.m[3][2] +
		matrix.m[0][1] * matrix.m[2][2] * matrix.m[3][0] +
		matrix.m[0][2] * matrix.m[2][0] * matrix.m[3][1] -
		matrix.m[0][2] * matrix.m[2][1] * matrix.m[3][0] -
		matrix.m[0][0] * matrix.m[2][2] * matrix.m[3][1] -
		matrix.m[0][1] * matrix.m[2][0] * matrix.m[3][2] ) / det;

	result.m[3][2] = ( matrix.m[0][0] * matrix.m[1][2] * matrix.m[3][1] +
		matrix.m[0][1] * matrix.m[1][0] * matrix.m[3][2] +
		matrix.m[0][2] * matrix.m[1][1] * matrix.m[3][0] -
		matrix.m[0][2] * matrix.m[1][0] * matrix.m[3][1] -
		matrix.m[0][0] * matrix.m[1][1] * matrix.m[3][2] -
		matrix.m[0][1] * matrix.m[1][2] * matrix.m[3][0] ) / det;

	result.m[3][3] = ( matrix.m[0][0] * matrix.m[1][1] * matrix.m[2][2] +
		matrix.m[0][1] * matrix.m[1][2] * matrix.m[2][0] +
		matrix.m[0][2] * matrix.m[1][0] * matrix.m[2][1] -
		matrix.m[0][2] * matrix.m[1][1] * matrix.m[2][0] -
		matrix.m[0][0] * matrix.m[1][2] * matrix.m[2][1] -
		matrix.m[0][1] * matrix.m[1][0] * matrix.m[2][2] ) / det;

	return result;
}

//5.転置行列
Matrix4x4 TransposeMatrix(const Matrix4x4& matrix) {
	Matrix4x4 result;
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			result.m[i][j] = matrix.m[j][i];
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
				result.m[i][j] = 1.0f;
			} else {
				result.m[i][j] = 0.0f;
			}
		}
	}
	return result;
}

#pragma endregion

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
			Novice::ScreenPrintf(x + column * kColumnWight, y + row * kRowHeight, "%6.02f", matrix.m[row][column]);
		}
	}

	Novice::ScreenPrintf(x + 250, y, "%s", label);
}

// 平行移動行列を作成する関数
Matrix4x4 MakeTranslateMatrix(const Vector3& translate) {
	Matrix4x4 result = IdentityMatrix(); // 単位行列を初期化

	// 平行移動成分をセット
	result.m[3][0] = translate.x;
	result.m[3][1] = translate.y;
	result.m[3][2] = translate.z;

	return result;
}

// 拡大縮小行列を作成する関数
Matrix4x4 MakeScaleMatrix(const Vector3& scale) {
	// 単位行列を初期化
	Matrix4x4 result = IdentityMatrix();

	// 拡大縮小成分をセット
	result.m[0][0] = scale.x;
	result.m[1][1] = scale.y;
	result.m[2][2] = scale.z;

	return result;
}

//X軸回転行列(yz)
Matrix4x4 MakeRotateXMatrix(float radian) {
	// 単位行列で初期化
	Matrix4x4 result = IdentityMatrix();

	//行列の計算
	result.m[1][1] = std::cos(radian);
	result.m[1][2] = std::sin(radian);
	result.m[2][1] = -std::sin(radian);
	result.m[2][2] = std::cos(radian);

	return result;
}

//Y軸回転行列(zx)
Matrix4x4 MakeRotateYMatrix(float radian) {
	// 単位行列で初期化
	Matrix4x4 result = IdentityMatrix();

	//行列の計算
	result.m[0][0] = std::cos(radian);
	result.m[0][2] = -std::sin(radian);
	result.m[2][0] = std::sin(radian);
	result.m[2][2] = std::cos(radian);

	return result;
}

//Z軸回転行列(xy)
Matrix4x4 MakeRotateZMatrix(float radian) {
	// 単位行列で初期化
	Matrix4x4 result = IdentityMatrix();

	//行列の計算
	result.m[0][0] = std::cos(radian);
	result.m[0][1] = std::sin(radian);
	result.m[1][0] = -std::sin(radian);
	result.m[1][1] = std::cos(radian);

	return result;
}

// 座標変換を行う関数
Vector3 Transform(const Vector3& vector, const Matrix4x4& matrix) {
	Vector3 result;

	// 行列とベクトルの乗算
	result.x = vector.x * matrix.m[0][0] + vector.y * matrix.m[1][0] + vector.z * matrix.m[2][0] + matrix.m[3][0];
	result.y = vector.x * matrix.m[0][1] + vector.y * matrix.m[1][1] + vector.z * matrix.m[2][1] + matrix.m[3][1];
	result.z = vector.x * matrix.m[0][2] + vector.y * matrix.m[1][2] + vector.z * matrix.m[2][2] + matrix.m[3][2];
	float w = vector.x * matrix.m[0][3] + vector.y * matrix.m[1][3] + vector.z * matrix.m[2][3] + matrix.m[3][3];
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
///レンダリングパイプライン
///

#pragma region レンダリングパイプライン


//正射影行列
Matrix4x4 MakeOrthographicMatrix(float left, float top, float right, float bottom, float nearClip, float farClip) {
	// 単位行列で初期化
	Matrix4x4 result = IdentityMatrix();

	// 正射影平面の範囲から正射影行列を構築する
	result.m[0][0] = 2.0f / ( right - left );
	result.m[1][1] = 2.0f / ( top - bottom );
	result.m[2][2] = -2.0f / ( farClip - nearClip );
	result.m[3][0] = -( right + left ) / ( right - left );
	result.m[3][1] = -( top + bottom ) / ( top - bottom );
	result.m[3][2] = -( farClip + nearClip ) / ( farClip - nearClip );

	return result;

}

//透視投影行列
Matrix4x4 MakePerspectiveFovMatrix(float fovY, float aspectRatio, float nearClip, float farClip) {
	// 単位行列で初期化
	Matrix4x4 result = IdentityMatrix();

	// tanの逆説
	float cot = tanf(fovY / 2.0f);

	// 射影変換行列の要素を計算する
	result.m[0][0] = 1.0f / ( aspectRatio * cot );
	result.m[1][1] = 1.0f / cot;
	result.m[2][2] = farClip / ( farClip - nearClip );
	result.m[2][3] = 1.0f;
	result.m[3][2] = -nearClip * farClip / ( farClip - nearClip );
	result.m[3][3] = 0.0f;

	return result;
}

//ビューポート行列
Matrix4x4 MakeViewportMatrix(float left, float top, float width, float height, float minDepth, float maxDepth) {
	// 単位行列で初期化
	Matrix4x4 result = IdentityMatrix();

	// ビューポート変換行列の要素を計算する
	result.m[0][0] = width / 2.0f;
	result.m[1][1] = -height / 2.0f;
	result.m[2][2] = maxDepth - minDepth;
	result.m[3][0] = left + width / 2.0f;
	result.m[3][1] = top + height / 2.0f;
	result.m[3][2] = minDepth;

	return result;
}
#pragma endregion

///
///クロス積
///
Vector3 Cross(const Vector3& v1, const Vector3& v2) {
	Vector3 result{};
	result.x = v1.y * v2.z - v1.z * v2.y;
	result.y = v1.z * v2.x - v1.x * v2.z;
	result.z = v1.x * v2.y - v1.y * v2.x;
	return result;
}

///
///Gridの表示
///
void DrawGrid(const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix) {
	const float kGridHalfWidth = 2.0f; // Gridの半分の幅
	const uint32_t kSubDivision = 12;  // 分割数
	const float kGridEvery = ( kGridHalfWidth * 2.0f ) / float(kSubDivision); // 1つ分の長さ

	// 奥から手前へ線を順々に引いていく
	for (uint32_t xIndex = 0; xIndex <= kSubDivision; ++xIndex) {
		Vector3 screenVertices[2];
		Vector3 kLocalVertices[2];

		// グリッドの頂点をワールド座標系で設定
		kLocalVertices[0] = Vector3(-kGridHalfWidth + xIndex * kGridEvery, 0.0f, kGridHalfWidth);
		kLocalVertices[1] = Vector3(-kGridHalfWidth + xIndex * kGridEvery, 0.0f, -kGridHalfWidth);

		// ワールド座標をスクリーン座標に変換
		for (uint32_t i = 0; i < 2; ++i) {
			Vector3 ndcVertex = Transform(kLocalVertices[i], viewProjectionMatrix);
			screenVertices[i] = Transform(ndcVertex, viewportMatrix);
		}

		if (xIndex == kSubDivision / 2) {
			// ラインを描画
			Novice::DrawLine(int(screenVertices[0].x), int(screenVertices[0].y), int(screenVertices[1].x), int(screenVertices[1].y), 0x555555FF);
		} else {
			// ラインを描画
			Novice::DrawLine(int(screenVertices[0].x), int(screenVertices[0].y), int(screenVertices[1].x), int(screenVertices[1].y), 0xAAAAAAFF);
		}
	}

	// 左から右へと同じように引いていく
	for (uint32_t zIndex = 0; zIndex <= kSubDivision; ++zIndex) {
		Vector3 screenVertices[2];
		Vector3 kLocalVertices[2];

		// グリッドの頂点をワールド座標系で設定
		kLocalVertices[0] = Vector3(kGridHalfWidth, 0.0f, -kGridHalfWidth + zIndex * kGridEvery);
		kLocalVertices[1] = Vector3(-kGridHalfWidth, 0.0f, -kGridHalfWidth + zIndex * kGridEvery);

		// ワールド座標をスクリーン座標に変換
		for (uint32_t i = 0; i < 2; ++i) {
			Vector3 ndcVertex = Transform(kLocalVertices[i], viewProjectionMatrix);
			screenVertices[i] = Transform(ndcVertex, viewportMatrix);
		}
		if (zIndex == kSubDivision / 2) {
			Novice::DrawLine(int(screenVertices[0].x), int(screenVertices[0].y), int(screenVertices[1].x), int(screenVertices[1].y), 0x555555FF);
		} else {
			// ラインを描画
			Novice::DrawLine(int(screenVertices[0].x), int(screenVertices[0].y), int(screenVertices[1].x), int(screenVertices[1].y), 0xAAAAAAFF);
		}
	}
}


///
///Sphereを描画する関数
///
void DrawSphere(const Sphere& sphere, const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix, uint32_t color) {
	const uint32_t kSubDivision = 20;  // 緯度と経度の分割数を増やす
	const float kLatEvery = float(M_PI) / float(kSubDivision);  // 緯度の刻み幅
	const float kLonEvery = 2.0f * float(M_PI) / float(kSubDivision);  // 経度の刻み幅

	// 経度方向に分割
	for (uint32_t latIndex = 0; latIndex <= kSubDivision; ++latIndex) {
		float lat = float(M_PI) / 2.0f - latIndex * kLatEvery;  // 現在の緯度
		float sinLat = sin(lat);
		float cosLat = cos(lat);

		// 次の緯度
		float nextLat = float(M_PI) / 2.0f - ( latIndex + 1 ) * kLatEvery;
		float sinNextLat = sin(nextLat);
		float cosNextLat = cos(nextLat);

		// 経度の方向に分割
		for (uint32_t lonIndex = 0; lonIndex < kSubDivision; ++lonIndex) {
			float lon = lonIndex * kLonEvery;  // 現在の経度
			float nextLon = ( lonIndex + 1 ) * kLonEvery;  // 次の経度

			// 現在の経度での頂点

			Vector3 a = AddVector3(MultiplyVector3(Vector3{ cos(lon) * cosLat, sin(lon) * cosLat, sinLat }, sphere.radius), sphere.center);
			Vector3 b = AddVector3(MultiplyVector3(Vector3{ cos(nextLon) * cosLat, sin(nextLon) * cosLat, sinLat }, sphere.radius), sphere.center);
			Vector3 c = AddVector3(MultiplyVector3(Vector3{ cos(lon) * cosNextLat, sin(lon) * cosNextLat, sinNextLat }, sphere.radius), sphere.center);
			Vector3 d = AddVector3(MultiplyVector3(Vector3{ cos(nextLon) * cosNextLat, sin(nextLon) * cosNextLat, sinNextLat }, sphere.radius), sphere.center);

			// ワールド座標をスクリーン座標に変換
			Vector3 screenA = Transform(Transform(a, viewProjectionMatrix), viewportMatrix);
			Vector3 screenB = Transform(Transform(b, viewProjectionMatrix), viewportMatrix);
			Vector3 screenC = Transform(Transform(c, viewProjectionMatrix), viewportMatrix);
			Vector3 screenD = Transform(Transform(d, viewProjectionMatrix), viewportMatrix);

			// ラインを描画
			Novice::DrawLine(int(screenA.x), int(screenA.y), int(screenB.x), int(screenB.y), color);
			Novice::DrawLine(int(screenA.x), int(screenA.y), int(screenC.x), int(screenC.y), color);
			Novice::DrawLine(int(screenB.x), int(screenB.y), int(screenD.x), int(screenD.y), color);
			Novice::DrawLine(int(screenC.x), int(screenC.y), int(screenD.x), int(screenD.y), color);
		}
	}
}


// Windowsアプリでのエントリーポイント(main関数)
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {

	const uint32_t kWindowWidth = 1280;
	const uint32_t kWindowHeight = 720;

	// ライブラリの初期化
	Novice::Initialize(kWindowTitle, 1280, 720);

	// キー入力結果を受け取る箱
	char keys[256] = { 0 };
	char preKeys[256] = { 0 };

	///-------------------------------
	///変数の宣言
	///-------------------------------
	//クロス積の確認用
	Vector3 v1{ 1.2f,-3.9f,2.5f };
	Vector3 v2{ 2.8f,0.4f,-1.3f };
	Vector3 cross = Cross(v1, v2);
	Vector3ScreenPrintf(0, 0, cross, "Cross");

	//世界
	Vector3 translate{ 0.0f,0.0f,0.0f };
	Vector3 rotate{ 0.0f,0.0f,0.0f };

	//カメラ行列
	Vector3 cameraTranslate{ 0.0f,1.9f,-6.49f };
	Vector3 cameraRotate{ 0.26f,0.0f,0.0f };

	//円
	Sphere sphere{ {0.0f,0.0f,0.0f},1.0f };


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

		// 各行列の計算
		Matrix4x4 worldMatrix = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, rotate, translate);
		Matrix4x4 cameraMatrix = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, cameraRotate, cameraTranslate);
		Matrix4x4 viewMatrix = InverseMatrix(cameraMatrix);
		Matrix4x4 projectionMatrix = MakePerspectiveFovMatrix(0.45f, float(kWindowWidth) / float(kWindowHeight), 0.1f, 100.0f);
		//
		Matrix4x4 worldViewProjectionMatrix = MultiplyMatrix(worldMatrix, MultiplyMatrix(viewMatrix, projectionMatrix));
		//
		Matrix4x4 viewportMatrix = MakeViewportMatrix(0, 0, float(kWindowWidth), float(kWindowHeight), 0.0f, 1.0f);
		//

		///
		/// ↑更新処理ここまで
		///

		///
		/// ↓描画処理ここから
		///

		//Gridの描画
		DrawGrid(worldViewProjectionMatrix, viewportMatrix);

		//円の描画
		DrawSphere(sphere,worldViewProjectionMatrix, viewportMatrix,0x777777FF);

		ImGui::Begin("Window");
		ImGui::DragFloat3("CameraTranslate", &cameraTranslate.x, 0.01f);
		ImGui::DragFloat3("CameraRotate", &cameraRotate.x, 0.01f);
		ImGui::DragFloat3("SphereCenter", &sphere.center.x, 0.01f);
		ImGui::DragFloat3("SphereRadius", &sphere.radius, 0.01f);
		ImGui::End();



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
