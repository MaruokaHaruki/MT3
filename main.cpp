#include <Novice.h>
#include <math.h>
#include<Vector3.h>

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
void Matrix4x4ScreenPrintf(int x, int y, const Matrix4x4& matrix, const char* label) {
	for (int row = 0; row < 4; ++row) {
		for (int column = 0; column < 4; ++column) {
			Novice::ScreenPrintf(x + column * kColumnWight, y + row * kRowHeight, "%6.02f", matrix.m4x4[row][column]);
		}
	}

	Novice::ScreenPrintf(x + 200, y, "%s", label);
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

// 平行移動行列を作成する関数
Matrix4x4 MakeTranslateMatrix(const Vector3& translate) {
	Matrix4x4 result = IdentityMatrix(); // 単位行列を初期化

	// 平行移動成分をセット
	result.m4x4[0][3] = translate.x;
	result.m4x4[1][3] = translate.y;
	result.m4x4[2][3] = translate.z;

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
	result.x = vector.x * matrix.m4x4[0][0] + vector.y * matrix.m4x4[0][1] + vector.z * matrix.m4x4[0][2] + matrix.m4x4[0][3];
	result.y = vector.x * matrix.m4x4[1][0] + vector.y * matrix.m4x4[1][1] + vector.z * matrix.m4x4[1][2] + matrix.m4x4[1][3];
	result.z = vector.x * matrix.m4x4[2][0] + vector.y * matrix.m4x4[2][1] + vector.z * matrix.m4x4[2][2] + matrix.m4x4[2][3];

	return result;
}


// Windowsアプリでのエントリーポイント(main関数)
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {

	// ライブラリの初期化
	Novice::Initialize(kWindowTitle, 1280, 720);

	// キー入力結果を受け取る箱
	char keys[256] = {0};
	char preKeys[256] = {0};

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
