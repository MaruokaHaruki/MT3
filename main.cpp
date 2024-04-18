#include <Novice.h>
#include<Vector3.h>
#include <math.h>

const char kWindowTitle[] = "LE2B_23_MT3_00_01";

///
///関数宣言
///

///3次元ベクトル表示
static void Vector3ScreenPrint(int x, int y, Vector3 vector3, const char* label) {
	Novice::ScreenPrintf(x, y, "%4.2f", vector3.x);
	Novice::ScreenPrintf(x + 50, y, "%4.2f", vector3.y);
	Novice::ScreenPrintf(x + 100, y, "%4.2f", vector3.z);
	Novice::ScreenPrintf(x + 150, y, "%s", label);
};

///加算
Vector3 Add(const Vector3& v1, const Vector3& v2) {
	Vector3 result;
	result.x = v1.x + v2.x;
	result.y = v1.y + v2.y;
	result.z = v1.z + v2.z;

	return result;
};

///減算
Vector3 Subtract(const Vector3& v1, const Vector3& v2) {
	Vector3 result;
	result.x = v1.x - v2.x;
	result.y = v1.y - v2.y;
	result.z = v1.z - v2.z;

	return result;
};

///スカラー倍
Vector3 Multiply(float scalar, const Vector3& v) {
	Vector3 result;
	result.x = scalar * v.x;
	result.y = scalar * v.y;
	result.z = scalar * v.z;

	return result;
};

///内積
float Dot(const Vector3& v1, const Vector3& v2) {
	float result;
	result = v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;

	return result;
};

///長さ(ノルム)
float Lenght(const Vector3& v) {
	float result;
	result = static_cast<float>(sqrt(v.x * v.x + v.y * v.y + v.z * v.z));

	return result;
};
///正規化
Vector3 Normalize(const Vector3& v) {
	float length = Lenght(v);
	// 長さが0の場合はそのまま返す（ゼロ除算を避ける）
	if (length == 0) {
		return v;
	}
	// 各成分を長さで割る
	return { v.x / length, v.y / length, v.z / length };
}



// Windowsアプリでのエントリーポイント(main関数)
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {

	// ライブラリの初期化
	Novice::Initialize(kWindowTitle, 1280, 720);

	// キー入力結果を受け取る箱
	char keys[256] = { 0 };
	char preKeys[256] = { 0 };

	///-----------------------------
	///変数の初期と代入
	///-----------------------------
	Vector3 v1{ 1.0f,3.0f,-5.0f };
	Vector3 v2{ 4.0f,-1.0f,2.0f };
	float k = { 4.0f };

	//計算結果の代入
	Vector3 resultAdd = Add(v1, v2);
	Vector3 resultSub = Subtract(v1, v2);
	Vector3 resultMult = Multiply(k, v1);
	float resultDot = Dot(v1, v2);
	float resultLenght = Lenght(v1);
	Vector3 resultNormalize = Normalize(v2);


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

		//加算
		Vector3ScreenPrint(10, 10, resultAdd, ":Add");
		//減産
		Vector3ScreenPrint(10, 30, resultSub, ":Subtract");
		//積
		Vector3ScreenPrint(10, 50, resultMult, ":Multiply");
		//内積
		Novice::ScreenPrintf(10, 70, "%4.2f :Dot", resultDot);
		//ノルム
		Novice::ScreenPrintf(10, 90, "%4.2f :Lenght", resultLenght);
		//正規化
		Vector3ScreenPrint(10, 110, resultNormalize, ":Normalize");

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
