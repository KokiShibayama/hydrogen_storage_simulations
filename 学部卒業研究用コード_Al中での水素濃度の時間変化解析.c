#include <stdio.h>
#include <math.h>

#define EV_TO_J 1.60219e-19 // eVをJに変換する係数

//E,θが入る別ファイルをコンパイル後に選べるようにコンパイル引数を設定
int main(int argc, char** argv) {

  
  //変数の設定
    //ns:単位面積当たりAl原子数, kb:ボルツマン定数, T:温度
    double n_s = 1.37e19; //大体の値。桁数がPdと同じ
    double k_B = 1.380649e-23; //(J/K)
    //double k_B = 8.617333262145e-5; //(eV/K)
    double T = 300; //単位:K
    //double p_H2 = 35e6; //35MPa(「論文：水素貯蔵」から引用)
    double p_H2 = 70e6; //70MPa(「論文：水素貯蔵」から引用)
    double m_H2 = 2.0e-3; //単位:kg
    double Γ = p_H2 / sqrt(2 * M_PI * m_H2 * k_B * T);
    //double Γ = 9702660584820180992.0;  //desmosで前もって計算
    //printf("n_s: %lf, k_B: %.17ef, T: %lf, p_H2: %lf, m_H2: %lf, Γ: %lf\n", n_s, k_B, T, p_H2, m_H2, Γ);  //確認用
  
  //E(エネルギー)の値を別ファイルから持ってくる

    FILE *fp1 = fopen(argv[1], "r");  //E記入ファイルを開く。名前は実行時に
    
    double E01, E02; //H_2↔sの活性化エネルギー
    double E1[7]; //Oi→Ti  //(値仮定)Oサイトが1~6,Tサイトが0~5で仮定している。
    double E2[7]; //Oi←Ti
    double E3[6]; //Tj→Oj+1
    double E4[6]; //Tj→Oj+1

    if (fp1 == NULL){
      printf("失敗(E)");  //ファイルが開けない（無い）場合の確認
    }
    
    fscanf(fp1, "%lf%lf", &E01, &E02);
    //printf("E01: %lf E02: %lf\n", E01, E02);  //確認用

    int i=1;
    int j=0;

          /*
          (使用しない)E[0]でT0↔O1を表現するように、TOとOTで分けたため。
          while (j<7 && fscanf(fp1, "%lf%lf%lf%lf", &E1[i], &E2[i], &E3[j], &E4[j]) == 4 ) { 
              printf("E1[%d]: %lf, E2[%d]: %lf, E3[%d]: %lf, E4[%d]: %lf", i, E1[i], i, E2[i], j, E3[j], j, E4[j]);  
          */

    while (i<7 && fscanf(fp1, "%lf%lf", &E1[i], &E2[i]) == 2 ){
        //printf("E1[%d]: %lf, E2[%d]: %lf\n", i, E1[i], i, E2[i]); //確認用
        i++;
    }
    while (j<6 && fscanf(fp1, "%lf%lf", &E3[j], &E4[j]) == 2 ){
        //printf("E3[%d]: %lf, E4[%d]: %lf\n", j, E3[j], j, E4[j]); //確認用
        j++;
    }
    fclose(fp1);
    //（実行メモ）Eが記入されているファイルを実行時にargv[1]として読み込めば全て代入される。Eの各順番には注意。
  
    //eVからJへの変換
    
      // E01, E02 の変換
      E01 *= EV_TO_J;
      E02 *= EV_TO_J;

      // 配列 E1[i], E2[i], E3[j], E4[j] の変換
      for (int i = 1; i < 7; i++) {
          E1[i] *= EV_TO_J;
          E2[i] *= EV_TO_J;
      }

      for (int j = 0; j < 6; j++) {
          E3[j] *= EV_TO_J;
          E4[j] *= EV_TO_J;
      }

      /*
      // 結果の出力(確認用)
      printf("E01 = %.10e J\n", E01);
      printf("E02 = %.10e J\n", E02);

      printf("E1: ");
      for (int i = 0; i < 6; i++) {
          printf("%.10e ", E1[i]);
      }
      printf("\n");

      printf("E2: ");
      for (int i = 0; i < 6; i++) {
          printf("%.10e ", E2[i]);
      }
      printf("\n");

      printf("E3: ");
      for (int j = 0; j < 6; j++) {
          printf("%.10e ", E3[j]);
      }
      printf("\n");

      printf("E4: ");
      for (int j = 0; j < 6; j++) {
          printf("%.10e ", E4[j]);
      }
      printf("\n");
      */

  //θ（初期濃度）の値を別ファイルから持ってくる

    FILE *fp2 = fopen(argv[2], "r");
    
    double θO[7]; //Oサイト層初期濃度([1]~[4])
    double θT[7]; //Tサイト層初期濃度([0]~[3])（注意θT[0]を表面サイトとして設定）

    if (fp2 == NULL){
      printf("失敗(θ)");  //ファイルが開けない（無い）場合の確認
    }

    //θOの代入
    i=1;
    j=0;
    
    while (i<7) {  //（注意）この条件はT,Oサイトの総数が同じだからできる。計算機戻って層の数分かったら見直す必要があるかもしれないから注意。
        fscanf(fp2, "%lf", &θT[j]);
        fscanf(fp2, "%lf", &θO[i]);
        //printf("θT[%d]: %lf, θO[%d]: %lf\n", j, θT[j], i, θO[i]);  //確認用
        i++;
        j++;
    }
    θT[j] = 0;
    //printf("θT[6]: %lf\n", θT[6]);
    //printf("θT[4]: %lf\n",  θT[4]);  //確認用
    fclose(fp2);
    //(実行メモ)θが記入されているファイルを実行時にargv[2]として読み込めば全て代入される。θは上(s=T0)から順に縦に並べる。

  //ν（ジャンプ頻度）の値を別ファイルから持ってくる

    FILE *fp3 = fopen(argv[3], "r");  //ν記入ファイルを開く。名前は実行時に
    
    double v0; //H2←s
    double v1[7]; //Oi→Ti
    double v2[7]; //Oi←Ti
    double v3[6]; //Ti→Oi+1
    double v4[6]; //Ti→Oi+1

    if (fp3 == NULL){
      //printf("失敗(v)");  //ファイルが開けない（無い）場合の確認
    }

    fscanf(fp3, "%lf", &v0);
    //printf("v0: %lf\n", v0);  //確認用

    i=1;
    j=0;

    while (i<7 && fscanf(fp3, "%lf%lf", &v1[i], &v2[i]) == 2 ){
        //printf("v1[%d]: %lf, v2[%d]: %lf\n", i, v1[i], i, v2[i]); //確認用
        i++;
    }
    while (j<6 && fscanf(fp3, "%lf%lf", &v3[j], &v4[j]) == 2 ){
        //printf("v3[%d]: %lf, v4[%d]: %lf\n", j, v3[j], j, v4[j]); //確認用
        j++;
    }

    fclose(fp3);
    //（実行メモ）vが記入されているファイルを実行時にargv[3]として読み込めば全て代入される。vの各順番には注意。

  //kの導出

    //(H2↔s(T0))の設定
      double k01, k02;
      //double E01;

      //k01=k(H2→s(T0)), k02=k(H2←s(T0))に相当
      k01 = (2.0 * Γ/n_s) * exp(-E01/(k_B * T));
      //k01 = exp(-E01/(k_B * T));
      k02 = v0 * exp(-E02/(k_B * T));
      //printf("k01: %lf, k02: %e\n", k01, k02);  //確認用.k02に関して、仮置きE02の値が大きいため計算結果がこれのみ非常に小さくなる。それが表示されるように%eにしている。

    //(Oi→Ti)の設定
      double k1[7];
      //double v1[6];
      //double E1[6];
      //a1[1]=a(O1→T1)に相当
      for (i=1; i<=6; i++){
        k1[i] = v1[i] * exp(-E1[i]/(k_B * T));
        //printf("k1[%d]: %lf\n", i, k1[i]);  //確認用
      }

    //(Oi←Ti)の設定
      double k2[7];
      //double v2[6];
      //double E2[6];
      //a1[1]=a1(O1→T1)に相当
      for (i=1; i<=6; i++){
        k2[i] = v2[i] * exp(-E2[i]/(k_B * T));
        //printf("k2[%d]: %lf\n", i, k2[i]);  //確認用
      }

    //(Ti→Oi+1)の設定
      double k3[6];
      //double v3[6];
      //double E3[6];
      //a3[0]=a3(s(T0)→O1), a3[1]=a3(T1→O2)に相当
      for (i=0; i<=5; i++){
        k3[i] = v3[i] * exp(-E3[i]/(k_B * T));
        //printf("k3[%d]: %lf\n", i, k3[i]);  //確認用
      }

    //(Ti←Oi+1)の設定
      double k4[6];
      //double v4[6];
      //double E4[6];
      //a4[0]=a3(s(T0)←O1), a4[1]=a3(T1←O2)に相当
      for (i=0; i<=5; i++){
        k4[i] = v4[i] * exp(-E4[i]/(k_B * T));
       // printf("k4[%d]: %lf\n", i, k4[i]);  //確認用
      }

  /*
  //標準出力からファイルに出力先を変えられるかの確認
  FILE *fp4 = freopen("output.txt", "w", stdout);  //出力を標準出力ではなくファイルに変更

  if (fp4 == NULL) {  //開けていない場合の確認
    printf("ファイルを開けませんでした\n");
    return 1;
    }
  
  printf("出力先変更成功");

  fflush(stdout);

  fclose(fp4);
  */

  //ループ部分
    
    //t(時間範囲の入力or定義)
    double t = 1.94e-12 ;  //(s)

    //試行回数の設定

    long n;

    FILE *fp4 = fopen("output.csv", "w");  //エクセルに読み込ませやすいようにcsvに
    //FILE *fp4 = freopen("output.csv", "w", stdout);  //出力を標準出力ではなくファイルに変更

    if (fp4 == NULL) {  //開けていない場合の確認
    printf("ファイルを開けませんでした\n");
    return 1;
    }
    
    /*
    if (fp4 == NULL) {
      perror("ファイルのオープンに失敗しました");
    }
    */

    /*
    // 標準出力をファイルにリダイレクト
    freopen("output.txt", "w", stdout);
    */

    fprintf(fp4, "n,θs\n");
    //printf("n,T[0],O[1],T[1],O[2],T[2],O[3],T[3],O[4],T[4],O[5],T[5],O[6]\n");

    //for (n=1; n<= 2 ; n++){  //検証
    for (n=1; n<= 10000000000000 ; n++){  //試行回数の設定  n=10^13

        //濃度の設定時間内での変化
        
        double Δθs;  //k01,k02使用のため表面のみ分けて考える。

        Δθs = t * ((k01*(1-θT[0])*(1-θT[0]) - k02*pow(θT[0], 2)) - (k3[0]*θT[0]*(1-θO[1]) - k4[0]*θO[1]*(1-θT[0])));  //使用変数:k(H2↔s), k(s↔O1), θT[0]

        double ΔθT[7], ΔθO[7];

        for (i=0; i<=5; i++){
          
          ΔθO[i+1] = t * ((k3[i]*θT[i]*(1-θO[i+1]) - k4[i]*θO[i+1]*(1-θT[i])) - (k1[i+1]*θO[i+1]*(1-θT[i+1]) - k2[i+1]*θT[i+1]*(1-θO[i+1])));
          
        }

        for (j=0; j<=5; j++){  //(注意)上の定義ではi=1ループでk[1]=O1関係というようにサイト番号と一致するようにしていた。ここではi=0ループにしてθs=θT0を使用するようにしているのでiの扱いに注意すること。

          ΔθT[j+1] = t * ((k1[j+1]*θO[j+1]*(1-θT[j+1]) - k2[j+1]*θT[j+1]*(1-θO[j+1])) - (k3[j+1]*θT[j+1]*(1-θO[j+2]) - k4[j+1]*θO[j+2]*(1-θT[j+1])));
                 
        }
        
        //最終層のみ次層を考慮しないため分ける。ここでは仮定でO層を採集にしているので、O6になる。
        
        //ΔθO[6] =  t * (k3[5]*θT[5]*(1-θO[6]) - k4[5]*θO[6]*(1-θT[5]));
        ΔθT[6] =  t * (k1[4]*θO[4]*(1-θT[4]) - k2[4]*θT[4]*(1-θO[4]));

        //変化分を足してt後の濃度の値導出       
        θT[0] = θT[0] + Δθs;
        //printf("θT[0]: %e\n", θT[0]);  //確認用

        for (i=1; i<=6; i++){
          θO[i] = θO[i] + ΔθO[i];
          //printf("θO[%d]: %e\n", i, θO[i]);  //確認用
        }

        for (j=1; j<=6; j++){
          θT[j] = θT[j] + ΔθT[j];
          //printf("θT[%d]: %e\n", j, θT[j]);  //確認用
        }
        
        //一定時間ごとに結果を記録
        
        if(n == 1 || (n % 10000000) == 0){  //10^7の倍数の時のみ
        //if(n==1 || n==1e1 || n==1e2 || n==1e3 || n==1e4 || n==1e5 || n==1e6 || n==1e7 || n==1e8 || n==1e9 || n==1e10 || n==1e11 || n==1e12 || n==1e13 || n==1e14){
        //if(n % 100 == 0){
            //fprintf(fp4, "%11d,%e,%e,%e,%e,%e,%e,%e,%e,%e\n", n, θT[0], θO[1], θT[1], θO[2], θT[2], θO[3], θT[3], θO[4], θT[4]);
            fprintf(fp4, "%ld,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n", n, θT[0], θO[1], θT[1], θO[2], θT[2], θO[3], θT[3], θO[4], θT[4], θO[5], θT[5], θO[6], θT[6]);
            fflush(fp4);    //強制終了対策

            printf("n=%ld, T0=%e,O1=%e,T1=%e,O2=%e,T2=%e,O3=%e,T3=%e,O4=%e,T4=%e,O5=%e,T5=%e,O6=%e,T6=%e\n", n, θT[0], θO[1], θT[1], θO[2], θT[2], θO[3], θT[3], θO[4], θT[4], θO[5], θT[5], θO[6], θT[6]);
            //printf("θT[0]:%e, θO[1]:%e\n", θT[0], θO[1]);
            printf("動いてるよ\n");
        }      
    }
      
    fclose(fp4);
    printf("データの書き出し完終わり\n");

  return 0;
}
