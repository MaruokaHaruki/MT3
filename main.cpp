#include <Novice.h>
#include<math.h>


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
//4x4行列の数値表示
static const int kRowHeight = 20;
static const int kColumnWight = 60;
void Matrix4x4ScreenPrintf(int x, int y, const Matrix4x4& matrix) {
	for (int row = 0; row < 4; ++row) {
		for (int column = 0; column < 4; ++column) {
			Novice::ScreenPrintf(x + column * kColumnWight, y + row * kRowHeight, "%6.02f", matrix.m4x4[row][column]);
		}
	}
}

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






// Windowsアプリでのエントリーポイント(main関数)
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {

	// ライブラリの初期化
	Novice::Initialize(kWindowTitle, 1280, 720);

	// キー入力結果を受け取る箱
	char keys[256] = { 0 };
	char preKeys[256] = { 0 };


	///-------------------------------
	///変数の宣言・初期化
	///-------------------------------
	Matrix4x4 m1 = { 3.2f,0.7f,9.6f,4.4f,
					5.5f,1.3f,7.8f,2.1f,
					6.9f,8.0f,2.6f,1.0f,
					0.5f,7.2f,5.1f,3.3f };

	Matrix4x4 m2 = { 4.1f,0.7f,9.6f,4.4f,
					8.8f,0.6f,9.9f,7.7f,
					1.1f,5.5f,6.6f,0.0f,
					3.3f,9.9f,8.8f,2.2f };


    Matrix4x4 resultAdd = AddMatrix(m1,m2);
    Matrix4x4 resultMult = MultiplyMatrix(m1, m2);
    Matrix4x4 resultSub = SubtractMatrix(m1, m2);
    Matrix4x4 inverseM1 = InverseMatrix(m1);
    Matrix4x4 inverseM2 = InverseMatrix(m2);
    Matrix4x4 transposeM1 = TransposeMatrix(m1);
    Matrix4x4 transposeM2 = TransposeMatrix(m2);
    Matrix4x4 identity = IdentityMatrix();


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
