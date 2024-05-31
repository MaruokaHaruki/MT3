#define NOMINMAX
#include <Novice.h>
#include <algorithm>
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

struct Line {
	Vector3 origin;
	Vector3 diff;
};

struct Ray {
	Vector3 origin;
	Vector3 diff;
};

struct Segment {
	Vector3 origin;
	Vector3 diff;
};

//平面
struct Plane {
	Vector3 normal;//法線
	float distance;//距離
};

struct Triangle {
	Vector3 vertics[3];
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

// ベクトルの大きさ（長さ）を計算する関数
float Magnitude(const Vector3& v) {
	return std::sqrt(Dot(v, v));
}

// 2つのベクトル間の距離を計算する関数
float Distance(const Vector3& a, const Vector3& b) {
	// ベクトルaとベクトルbの各成分の差を計算し、その平方を求める
	// 各成分の平方を合計し、その平方根を取ることで距離を計算する
	return std::sqrt(( a.x - b.x ) * ( a.x - b.x ) + ( a.y - b.y ) * ( a.y - b.y ) + ( a.z - b.z ) * ( a.z - b.z ));
}

// ベクトルを正規化
Vector3 Normalize(const Vector3& v) {
	float mag = Magnitude(v);
	if (mag != 0.0f) {
		return { v.x / mag, v.y / mag, v.z / mag };
	}
	// ゼロベクトルの場合はそのまま返す
	return v;
}

// Length 関数の実装
float Length(const Vector3& v) {
	return std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
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

// 行列の余因子行列を計算するヘルパー関数
Matrix4x4 CofactorMatrix(const Matrix4x4& matrix);

// 余因子行列の計算に必要なサブ行列を計算するヘルパー関数
float Minor(const Matrix4x4& matrix, int row, int col);

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

	// 行列式が0の場合は逆行列は存在しない


	Matrix4x4 cofactorMatrix = CofactorMatrix(matrix);
	Matrix4x4 adjugateMatrix;
	// 余因子行列の転置を求める
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			adjugateMatrix.m[i][j] = cofactorMatrix.m[j][i];
		}
	}

	Matrix4x4 inverseMatrix;
	// 逆行列を求める
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			inverseMatrix.m[i][j] = adjugateMatrix.m[i][j] / det;
		}
	}

	return inverseMatrix;
}

Matrix4x4 CofactorMatrix(const Matrix4x4& matrix) {
	Matrix4x4 cofactorMatrix;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			float minor = Minor(matrix, i, j);
			// 余因子行列の符号付き小行列行列式
			cofactorMatrix.m[i][j] = ( ( i + j ) % 2 == 0 ? 1 : -1 ) * minor;
		}
	}
	return cofactorMatrix;
}

float Minor(const Matrix4x4& matrix, int row, int col) {
	float subMatrix[3][3];
	int subRow = 0;
	for (int i = 0; i < 4; i++) {
		if (i == row) continue;
		int subCol = 0;
		for (int j = 0; j < 4; j++) {
			if (j == col) continue;
			subMatrix[subRow][subCol] = matrix.m[i][j];
			subCol++;
		}
		subRow++;
	}

	return subMatrix[0][0] * ( subMatrix[1][1] * subMatrix[2][2] - subMatrix[1][2] * subMatrix[2][1] ) -
		subMatrix[0][1] * ( subMatrix[1][0] * subMatrix[2][2] - subMatrix[1][2] * subMatrix[2][0] ) +
		subMatrix[0][2] * ( subMatrix[1][0] * subMatrix[2][1] - subMatrix[1][1] * subMatrix[2][0] );
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

Matrix4x4 MakeRotateMatrix(const Vector3& rotate) {
	Matrix4x4 rx = MakeRotateXMatrix(rotate.x);
	Matrix4x4 ry = MakeRotateYMatrix(rotate.y);
	Matrix4x4 rz = MakeRotateZMatrix(rotate.z);

	// 回転行列の掛け合わせ順序はY, X, Zとする
	return MultiplyMatrix(ry, MultiplyMatrix(rx, rz));
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
/// 
#pragma region Gridの表示

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
#pragma endregion

///
///Sphereを描画する関数
///

void DrawSphere(const Sphere& sphere, const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix, uint32_t color) {
	const uint32_t kSubDivision = 16;
	const float kLatEvery = float(M_PI) / float(kSubDivision);
	const float kLonEvery = 2.0f * float(M_PI) / float(kSubDivision);

	for (uint32_t latIndex = 0; latIndex <= kSubDivision; ++latIndex) {
		float lat = float(M_PI) / 2.0f - latIndex * kLatEvery;
		float sinLat = sin(lat);
		float cosLat = cos(lat);

		float nextLat = float(M_PI) / 2.0f - ( latIndex + 1 ) * kLatEvery;
		float sinNextLat = sin(nextLat);
		float cosNextLat = cos(nextLat);

		for (uint32_t lonIndex = 0; lonIndex < kSubDivision; ++lonIndex) {
			float lon = lonIndex * kLonEvery;
			float nextLon = ( lonIndex + 1 ) * kLonEvery;

			Vector3 a = { cos(lon) * cosLat, sin(lon) * cosLat, sinLat };
			Vector3 b = { cos(nextLon) * cosLat, sin(nextLon) * cosLat, sinLat };
			Vector3 c = { cos(lon) * cosNextLat, sin(lon) * cosNextLat, sinNextLat };
			Vector3 d = { cos(nextLon) * cosNextLat, sin(nextLon) * cosNextLat, sinNextLat };

			a = AddVector3(MultiplyVector3(a, sphere.radius), sphere.center);
			b = AddVector3(MultiplyVector3(b, sphere.radius), sphere.center);
			c = AddVector3(MultiplyVector3(c, sphere.radius), sphere.center);
			d = AddVector3(MultiplyVector3(d, sphere.radius), sphere.center);

			Vector3 screenA = Transform(Transform(a, viewProjectionMatrix), viewportMatrix);
			Vector3 screenB = Transform(Transform(b, viewProjectionMatrix), viewportMatrix);
			Vector3 screenC = Transform(Transform(c, viewProjectionMatrix), viewportMatrix);
			Vector3 screenD = Transform(Transform(d, viewProjectionMatrix), viewportMatrix);

			// ここでラインを描画する関数を呼び出す
			// ラインを描画
			Novice::DrawLine(int(screenA.x), int(screenA.y), int(screenB.x), int(screenB.y), color);
			Novice::DrawLine(int(screenA.x), int(screenA.y), int(screenC.x), int(screenC.y), color);
			Novice::DrawLine(int(screenB.x), int(screenB.y), int(screenD.x), int(screenD.y), color);
			Novice::DrawLine(int(screenC.x), int(screenC.y), int(screenD.x), int(screenD.y), color);
		}
	}
}


///
///正射影ベクトル
/// 
Vector3 Project(const Vector3& v1, const Vector3& v2) {
	float dotProduct = Dot(v1, v2);
	float magnitudeSquared = Magnitude(v2) * Magnitude(v2);
	return MultiplyVector3(v2, dotProduct / magnitudeSquared);
}

///
///最近接点
///
Vector3 ClossPoint(const Vector3& point, const Segment& segment) {
	// 線分の方向ベクトルを計算
	Vector3 segmentDirection = SubtractVector3(segment.diff, segment.origin);

	// 点から線分の始点までのベクトルを計算
	Vector3 pointToOrigin = SubtractVector3(point, segment.origin);

	// 点から始点までのベクトルを線分の方向ベクトルに射影
	Vector3 projection = Project(pointToOrigin, segmentDirection);

	// 射影ベクトルと線分方向ベクトルの内積を線分方向ベクトルの長さの2乗で割ることでtを求める
	float t = Dot(projection, segmentDirection) / Dot(segmentDirection, segmentDirection);

	// tを0と1の間に制限（線分の範囲外の場合に線分の端点にクランプ）
	t = std::max(0.0f, std::min(t, 1.0f));

	// 線分方向ベクトルにtを掛けたベクトルを返す
	return AddVector3(segment.origin, MultiplyVector3(segmentDirection, t));
}

///
///平面の描画
///
Vector3 Perpendicular(const Vector3& vector) {
	if (vector.x != 0.0f || vector.y != 0.0f) {
		return { -vector.y , vector.x,0.0f };
	}
	return { 0.0f,-vector.z,vector.y };
}

void DrawPlane(const Plane& plane, const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix, uint32_t color) {
	//1
	Vector3 center = MultiplyVector3(plane.normal, plane.distance);
	Vector3 perpendiculars[4];
	//2
	perpendiculars[0] = Normalize(Perpendicular(plane.normal));
	//3
	perpendiculars[1] = { -perpendiculars[0].x,-perpendiculars[0].z };
	//4
	perpendiculars[2] = Cross(plane.normal, perpendiculars[0]);
	//5
	perpendiculars[3] = { -perpendiculars[2].x,-perpendiculars[2].y,-perpendiculars[2].z };
	//6
	Vector3 points[4];
	for (int32_t index = 0; index < 4; ++index) {
		Vector3 extend = MultiplyVector3(perpendiculars[index], 2.0f);
		Vector3 point = AddVector3(center, extend);
		points[index] = Transform(Transform(point, viewProjectionMatrix), viewportMatrix);
	}
	//描画
	Novice::DrawLine(static_cast<int>( points[0].x ), static_cast<int>( points[0].y ), static_cast<int>( points[1].x ), static_cast<int>( points[1].y ), color);
	Novice::DrawLine(static_cast<int>( points[0].x ), static_cast<int>( points[0].y ), static_cast<int>( points[2].x ), static_cast<int>( points[2].y ), color);
	Novice::DrawLine(static_cast<int>( points[0].x ), static_cast<int>( points[0].y ), static_cast<int>( points[3].x ), static_cast<int>( points[3].y ), color);

	Novice::DrawLine(static_cast<int>( points[1].x ), static_cast<int>( points[1].y ), static_cast<int>( points[2].x ), static_cast<int>( points[2].y ), color);
	Novice::DrawLine(static_cast<int>( points[1].x ), static_cast<int>( points[1].y ), static_cast<int>( points[3].x ), static_cast<int>( points[3].y ), color);

	Novice::DrawLine(static_cast<int>( points[2].x ), static_cast<int>( points[2].y ), static_cast<int>( points[3].x ), static_cast<int>( points[3].y ), color);
	Novice::DrawLine(static_cast<int>( points[3].x ), static_cast<int>( points[3].y ), static_cast<int>( points[0].x ), static_cast<int>( points[0].y ), color);

}


///
///三角形の描画
///
void DrawTriangle(const Triangle& triangle, const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix, uint32_t color) {
	// 三角形の各頂点を変換する
	Vector3 transformedVertices[3];
	for (int32_t index = 0; index < 3; ++index) {
		transformedVertices[index] = Transform(Transform(triangle.vertics[index], viewProjectionMatrix), viewportMatrix);
	}

	// 描画
	Novice::DrawLine(static_cast<int>( transformedVertices[0].x ), static_cast<int>( transformedVertices[0].y ), static_cast<int>( transformedVertices[1].x ), static_cast<int>( transformedVertices[1].y ), color);
	Novice::DrawLine(static_cast<int>( transformedVertices[0].x ), static_cast<int>( transformedVertices[0].y ), static_cast<int>( transformedVertices[2].x ), static_cast<int>( transformedVertices[2].y ), color);
	Novice::DrawLine(static_cast<int>( transformedVertices[1].x ), static_cast<int>( transformedVertices[1].y ), static_cast<int>( transformedVertices[2].x ), static_cast<int>( transformedVertices[2].y ), color);
}



///
///球体の衝突判定
///
bool IsCollision(const Sphere& s1, const Sphere& s2) {
	// 2つの球体の中心間の距離を計算
	float centerDistance = Distance(s1.center, s2.center);

	// 2つの球体の半径の合計
	float radiusSum = s1.radius + s2.radius;

	// 中心間の距離が半径の合計以下であれば衝突していると判定
	return centerDistance <= radiusSum;
}

///
///球体と平面の衝突判定
/// 
bool IsSphere2PlaneCollision(const Sphere& sphere, const Plane& plane) {
	// 球体の中心から平面までの距離を計算
	float distance = Dot(plane.normal, sphere.center) - plane.distance;
	// 距離が球体の半径以内であれば衝突している
	return fabs(distance) <= sphere.radius;
}

///
///線と平面の衝突判定
///

bool IsLine2Sphere(const Segment& segment, const Plane& plane) {
	//垂直判定を行うために、法線と線の内積を求める
	float dot = Dot(plane.normal, segment.diff);

	//垂直=並行であるので、衝突しているはずがない
	if (dot == 0.0f) {
		return false;
	}

	//tを求める
	float t = ( plane.distance - Dot(segment.origin, plane.normal) ) / dot;

	// tの値と線の種類によって衝突しているかを判定する
	if (t >= 0.0f && t <= 1.0f) {
		// 線分が平面と交差している
		return true;
	} else {
		// 線分が平面と交差していない
		return false;
	}
}

///
///三角形と線の衝突判定
///
bool IsTriangle2Line(const Triangle& triangle, const Segment& segment) {
	Vector3 edge1 = triangle.vertics[1] - triangle.vertics[0];
	Vector3 edge2 = triangle.vertics[2] - triangle.vertics[0];
	Vector3 h = Cross(segment.diff, edge2);
	float a = Dot(edge1, h);
	if (fabs(a) < 1e-5) {
		// 線分は三角形の平面と平行
		return false;
	}

	float f = 1.0f / a;
	Vector3 s = segment.origin - triangle.vertics[0];
	float u = f * Dot(s, h);
	if (u < 0.0f || u > 1.0f) {
		return false;
	}

	Vector3 q = Cross(s, edge1);
	float v = f * Dot(segment.diff, q);
	if (v < 0.0f || u + v > 1.0f) {
		return false;
	}

	// tは線分のパラメータ、交差点の位置を計算
	float t = f * Dot(edge2, q);
	if (t < 0.0f || t > 1.0f) {
		// 線分の範囲外に交差点がある
		return false;
	}

	// 交差点は三角形の内部にある
	return true;
}

///
///カメラの位置
///
Matrix4x4 LookAt(const Vector3& eye, const Vector3& target, const Vector3& up) {
	Vector3 zaxis = Normalize(target - eye);    // 前方向ベクトル
	Vector3 xaxis = Normalize(Cross(up, zaxis)); // 右方向ベクトル
	Vector3 yaxis = Cross(zaxis, xaxis);        // 上方向ベクトル

	Matrix4x4 viewMatrix = {
		xaxis.x, yaxis.x, zaxis.x, 0,
		xaxis.y, yaxis.y, zaxis.y, 0,
		xaxis.z, yaxis.z, zaxis.z, 0,
		-Dot(xaxis, eye), -Dot(yaxis, eye), -Dot(zaxis, eye), 1
	};

	return viewMatrix;
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

	// カメラ行列
	Vector3 cameraTranslate{ 0.0f, 1.9f, -10.49f };
	Vector3 cameraRotate{ 0.26f, 0.0f, 0.0f };
	Sphere cameraTarget;
	cameraTarget.center = { 0.0f, 0.0f, 0.0f }; // カメラのターゲットポイント
	cameraTarget.radius = 0.01f;
	
	int lastMouseX = 0;
	int lastMouseY = 0;
	int mouseX = 0;
	int mouseY = 0;
	bool IsDebugCameraActive = false;

	//円
	Sphere sphere1{ {0.0f,0.0f,0.0f},1.0f };
	Sphere sphere2{ {1.0f,1.0f,0.0f},1.0f };

	//平面
	Plane plane{ {0.0f,1.0f,0.0f}, { 0.0f } };


	//点
	Segment segment{ {-2.0f,-1.0f,0.0f},{3.0f,2.0f,2.0f} };
	//Vector3 point{ -1.5f,0.6f,0.6f };

	//正射影ベクトルの計算
	//Vector3 project = Project(SubtractVector3(point, segment.origin), segment.diff);
	////最近接点
	//Vector3 clossPoint = ClossPoint(point, segment);

	//Sphere pointSphere{ point,0.01f };
	//Sphere clossPointSphere{ clossPoint,0.01f };


	//三角形
	Triangle triangle;
	triangle.vertics[0] = { 0.0f,2.0f,0.0f };
	triangle.vertics[1] = { 2.0f,-2.0f,0.0f };
	triangle.vertics[2] = { -2.0f,-2.0f,0.0f };

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

		/// ===デバックカメラ起動=== ///
		if (Novice::CheckHitKey(DIK_SPACE)) {
			if (IsDebugCameraActive) {
				IsDebugCameraActive = false;
			} else {
				IsDebugCameraActive = true;
			}
		}


		if (IsDebugCameraActive) {
			Novice::GetMousePosition(&mouseX, &mouseY);

			if (Novice::IsPressMouse(0)) {
				

				// マウスの移動量を計算
				int deltaX = mouseX - lastMouseX;
				int deltaY = mouseY - lastMouseY;

				// カメラの回転を更新
				float rotationSpeed = 0.005f;
				cameraRotate.y += deltaX * rotationSpeed;
				cameraRotate.x += deltaY * rotationSpeed;

				// カメラの位置をターゲットポイントの周りに回転
				float distance = Length(cameraTranslate - cameraTarget.center);
				Matrix4x4 rotationMatrix = MakeRotateMatrix(cameraRotate);
				Vector3 offset = { 0.0f, 0.0f, -distance };
				cameraTranslate = cameraTarget.center + Transform(offset, rotationMatrix);


			}
			// マウスの位置を更新
			lastMouseX = mouseX;
			lastMouseY = mouseY;

		}



		// 各行列の計算
		Matrix4x4 worldMatrix = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, rotate, translate);
		Matrix4x4 cameraMatrix = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, cameraRotate, cameraTranslate);


		Matrix4x4 viewMatrix = LookAt(cameraTranslate, cameraTarget.center, { 0.0f, 1.0f, 0.0f });
		Matrix4x4 projectionMatrix = MakePerspectiveFovMatrix(0.45f, float(kWindowWidth) / float(kWindowHeight), 0.1f, 100.0f);
		Matrix4x4 worldViewProjectionMatrix = MultiplyMatrix(worldMatrix, MultiplyMatrix(viewMatrix, projectionMatrix));
		Matrix4x4 viewportMatrix = MakeViewportMatrix(0, 0, float(kWindowWidth), float(kWindowHeight), 0.0f, 1.0f);

		//Matrix4x4 viewMatrix = InverseMatrix(cameraMatrix);
		//Matrix4x4 projectionMatrix = MakePerspectiveFovMatrix(0.45f, float(kWindowWidth) / float(kWindowHeight), 0.1f, 100.0f);
		////ワールドビューマトリックス
		//Matrix4x4 worldViewProjectionMatrix = MultiplyMatrix(worldMatrix, MultiplyMatrix(viewMatrix, projectionMatrix));
		////ビューポイント
		//Matrix4x4 viewportMatrix = MakeViewportMatrix(0, 0, float(kWindowWidth), float(kWindowHeight), 0.0f, 1.0f);
		//



		ImGui::Begin("Window");

		// カメラ関連
		ImGui::Text("Camera Settings");
		ImGui::Text("Press SPACE is DebugCameraActive");
		ImGui::Separator();
		ImGui::Spacing();
		ImGui::DragFloat3("Camera Translate", &cameraTranslate.x, 0.01f);
		ImGui::DragFloat3("Camera Rotate", &cameraRotate.x, 0.01f);
		ImGui::DragFloat3("Camera Target", &cameraTarget.center.x, 0.01f); // ターゲットポイントの調整

		// オブジェクト関連
		ImGui::Spacing();
		ImGui::Text("Object Settings");
		ImGui::Separator();
		ImGui::Spacing();
		// 球体1
		ImGui::DragFloat3("Sphere 1 Center", &sphere1.center.x, 0.01f);
		ImGui::DragFloat("Sphere 1 Radius", &sphere1.radius, 0.01f);
		// 球体2
		ImGui::DragFloat3("Sphere 2 Center", &sphere2.center.x, 0.01f);
		ImGui::DragFloat("Sphere 2 Radius", &sphere2.radius, 0.01f);
		// 平面
		ImGui::DragFloat3("Plane Normal", &plane.normal.x, 0.01f);
		// NOTE: 法線を編集したらNormalizeをかけること。平面法線が単位ベクトル前提でアルゴリズムが組まれているため
		plane.normal = Normalize(plane.normal);
		ImGui::DragFloat("Plane Distance", &plane.distance, 0.01f);
		// 線
		ImGui::DragFloat3("Segment Origin", &segment.origin.x, 0.01f);
		ImGui::DragFloat3("Segment Diff", &segment.diff.x, 0.01f);
		// 三角形
		ImGui::DragFloat3("Triangle Vertex 1", &triangle.vertics[0].x, 0.01f);
		ImGui::DragFloat3("Triangle Vertex 2", &triangle.vertics[1].x, 0.01f);
		ImGui::DragFloat3("Triangle Vertex 3", &triangle.vertics[2].x, 0.01f);

		ImGui::End();

		///
		/// ↑更新処理ここまで
		///

		///
		/// ↓描画処理ここから
		///

		//平面の描画
		DrawPlane(plane, worldViewProjectionMatrix, viewportMatrix, WHITE);

		//中心の描画
		DrawSphere(cameraTarget, worldViewProjectionMatrix, viewportMatrix, RED);

		////円の描画
		//if (IsCollision(sphere1, sphere2)) {
		//	DrawSphere(sphere1, worldViewProjectionMatrix, viewportMatrix, RED);
		//} else {
		//	DrawSphere(sphere1, worldViewProjectionMatrix, viewportMatrix, WHITE);
		//}
		//DrawSphere(sphere2, worldViewProjectionMatrix, viewportMatrix, WHITE);

		//Gridの描画

		////塩と平面の接触判定
		//if (IsSphere2PlaneCollision(sphere1, plane)) {
		//	DrawSphere(sphere1, worldViewProjectionMatrix, viewportMatrix, RED);
		//} else {
		//	DrawSphere(sphere1, worldViewProjectionMatrix, viewportMatrix, WHITE);
		//}

		////三角形と線の判定
		////if (IsTriangle2Line(triangle, segment)) {
		////	三角形の描画
		////	DrawTriangle(triangle, worldViewProjectionMatrix, viewportMatrix, RED);
		////} else {
		////	DrawTriangle(triangle, worldViewProjectionMatrix, viewportMatrix, WHITE);
		////}
		DrawTriangle(triangle, worldViewProjectionMatrix, viewportMatrix, WHITE);

		DrawGrid(worldViewProjectionMatrix, viewportMatrix);


		//線分の描画
		Vector3 start = Transform(Transform(segment.origin, worldViewProjectionMatrix), viewportMatrix);
		Vector3 end = Transform(Transform(AddVector3(segment.origin, segment.diff), worldViewProjectionMatrix), viewportMatrix);

		////塩と平面の接触判定
		//if (IsLine2Sphere(segment, plane)) {
		//	Novice::DrawLine(int(start.x), int(start.y), int(end.x), int(end.y), RED);
		//} else {
		//	Novice::DrawLine(int(start.x), int(start.y), int(end.x), int(end.y), WHITE);
		//}

		//塩と平面の接触判定
		if (IsTriangle2Line(triangle, segment)) {
			Novice::DrawLine(int(start.x), int(start.y), int(end.x), int(end.y), RED);
		} else {
			Novice::DrawLine(int(start.x), int(start.y), int(end.x), int(end.y), WHITE);
		}


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
