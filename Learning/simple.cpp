#include <iostream>
#include <typeinfo>
#include <vector>
#include <random>
#include <fstream>
#include <sstream>

using namespace std;

int mult2(int x){
    return x*2;
}

vector<int> mult2(vector<int> x, int size){
    vector<int> result = {};
    for(auto elem: x){
        result.push_back(mult2(elem));
    }
    return result;
}

auto mult2(auto numbers){
    vector<int> result = {};
    for(auto elem: numbers){
        result.push_back(mult2(elem));
    }    
    return result;
}
int old(){
    char s[] = "42lasdgkn";

    //int i{std::atoi(s)};
    int i = atoi(s);
    //std::cout << i << endl;


    
    int arr2[4] = {1,5,2,3};

    double result[4] = {};
    
    auto vec = {1,2,3,4};
    vector<int> lolo = mult2(vec);

    // for(int i=0; i<lolo.size(); i++){
    //     cout << lolo[i] << endl;
    // }

    int arr1[5] = {4,2,7,12,10};
    int size = sizeof(arr1)/sizeof(arr1[0]);

    cout << size << endl;
    return 0;
}
auto compareA(int x, int y){
    return x>y;
}
auto compareAll(){
    int xs[2] = {1,2};
    int ys[2] = {0,20};

    bool results[2] = {};
    for(int i=0; i<2; i++){
        results[i] = compareA(xs[i],ys[i]);
    }

    for(auto r: results){
        cout << r << endl;
    }
    return *xs;
}


auto def(int size){
    double *ns = new double[size];
    for(int i=0; i<size; i++){
        ns[i] = i;
    }
    return ns;
}

double compare(double x, double y){
    return x>y;
}
auto compare(auto xs, auto ys, int size){
    int *results = new int[size];
    for(int i=0; i<size; i++){
        if(xs[i] > ys[i]){
            results[i] = 1;
        }
        else{
            results[i] = 0;
        }
    }
    return results;
}

double awesome_sig(double z){
    double z_4 = z/4;
    cout << "x/4   " <<z_4 << endl;

    double term1 = z_4+0.5;
    cout << "x/4 + 0.5    " <<term1 << endl;

    //double term2 = min(1,term1);
    //cout << "min(1,term1)   " <<term2 << endl;

    //double sig = max(0,term2);
    auto sig=0;
    return sig;
}

double** Transpose(double **matrix, int n, int m){
    double** transpose = new double*[m];
    for(int i=0; i<m; i++){
        transpose[i] = new double[n]; //memory allocation baby!!
    }

    for(int i=0; i<n; ++i){
        for(int j=0; j<m; ++j){
            transpose[j][i] = matrix[i][j];
        }
    }
    return transpose;
}

int* multiplyMatrixVector(int** matrix, int* vec, int r, int n) {
    // Allocate memory for the result vector
    int* result = new int[r]();

    // Perform matrix-vector multiplication
    for (int i = 0; i < r; ++i) {
        for (int j = 0; j < n; ++j) {
            result[i] += matrix[i][j] * vec[j];
        }
    }

    return result;
}
double** multiplyMatrices(double** matrixA, double** matrixB, int nA, int mA, int mB) {
    // Allocate memory for the result matrix
    double** result = new double*[nA];
    for (int i = 0; i < nA; ++i) {
        result[i] = new double[mB]();
    }
    // Perform matrix multiplication
    for (int i = 0; i < nA; ++i) {
        for (int j = 0; j < mB; ++j) {
            for (int k = 0; k < mA; ++k) {
                result[i][j] += matrixA[i][k] * matrixB[k][j];
            }
        }
    }
    return result;
}

void increment(int *theta, int size){
    for(int i=0; i<size; i++){
        theta[i] += 1;
    }
}
int oldmain(){
    //cout << awesome_sig(0.75) << endl;

    random_device rd{};
    mt19937 gen{rd()};
 
    // values near the mean are the most likely
    // standard deviation affects the dispersion of generated values from the mean
    normal_distribution d{.0, 0.1}; //mean, std_dev
 
    // draw a sample from the normal distribution and round it to an integer
    auto rand_number = [&d, &gen] {return d(gen);};
    cout << d(gen) << endl;
    //auto random_int = [&d, &gen]{ return std::round(d(gen)); };

    
    vector<int> lololo = {1,2,3};
    for(int i=0; i<lololo.size(); ++i){
        cout << lololo[i] << endl;
    }

    int size = 3;
    int *theta = new int[size];




    int **matris = new int*[3];




    for(int i=0; i<size; i++){
        theta[i] = 10;
    }
    increment(theta, size);

    for(int i=0; i<size; i++){
        cout << theta[i] << " ";
    }
    cout <<endl;
    

    //cout << product[0] << endl;
    //cout << product[1] << endl;
    //matris = T(matris, n, m);
    // for(int i=0; i<n; ++i){
    //     for(int j=0; j<m; ++j){
    //         cout << matris[i][j] << " ";
    //     }
    //     cout << endl;
    // }
    return 0;
}



int read(){
    std::ifstream infile; 
    infile.open("/Users/dulumrae/Desktop/preprocessed_data.txt");   //burayı senin değiştirmen lazım!!!!!!!!!

    if (!infile.is_open()) {
        std::cerr << "Failed to open the file." << std::endl;
        
        std::cerr << "Error code: " << errno << " - " << std::strerror(errno) << std::endl; 

        return 1;
    }

    const int numRows = 30162; // number of rows
    const int numCols = 15;  // number of columns (target column dahil)

    double** X = new double*[numRows]; // Matrix X
    for (int i = 0; i < numRows; ++i) {
        X[i] = new double[numCols - 1]; 
    }
    double* y = new double[numRows]; // vector y
    std::string line;
    int row = 0;

    while (std::getline(infile, line) && row < numRows) {
        std::stringstream ss(line);
        std::string token;
        int col = 0;

        
        while (std::getline(ss, token, ',') && col < numCols - 1) { // , ile ayrılıyorlar
            X[row][col] = std::stod(token); // Convert string to double
            col++;
        }

        // target için bu kısım (sondaki sütun)
        if (col == numCols - 1 && std::getline(ss, token, ',')) {
            y[row] = std::stod(token);
        }

        row++;
    }

    infile.close();

    // Debnemek için yazdım
    std::cout << "Features (X):" << std::endl;
    for (int i = 0; i < numRows; ++i) {
        for (int j = 0; j < numCols - 1; ++j) {
            std::cout << X[i][j] << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "Target values (y):" << std::endl;
    for (int i = 0; i < numRows; ++i) {
        std::cout << y[i] << " ";
    }
    std::cout << std::endl;

    // Clean up
    for (int i = 0; i < numRows; ++i) {
        delete[] X[i];
    }
    delete[] X;
    delete[] y;

    return 0;
}


auto read_and_return(){
    const int numRows = 30162; // number of rows
    const int numCols = 15;  // number of columns (target column dahil)
    double** X_read = new double*[numRows]; // Matrix X
    for (int i = 0; i < numRows; ++i)   {X_read[i] = new double[numCols - 1]; }
    double* y_read = new double[numRows]; // vector y


    std::ifstream infile; 
    infile.open("/Users/dulumrae/Desktop/preprocessed_data.txt");   //burayı senin değiştirmen lazım!!!!!!!!!

    if (!infile.is_open()) {
        std::cerr << "Failed to open the file." << std::endl;
        
        std::cerr << "Error code: " << errno << " - " << std::strerror(errno) << std::endl; 

        return 1;
    }
    
    std::string line;
    int row = 0;
    while (std::getline(infile, line) && row < numRows) {
        std::stringstream ss(line);
        std::string token;
        int col = 0;

        
        while (std::getline(ss, token, ',') && col < numCols - 1) { // , ile ayrılıyorlar
            X_read[row][col] = std::stod(token); // Convert string to double
            col++;
        }

        // target için bu kısım (sondaki sütun)
        if (col == numCols - 1 && std::getline(ss, token, ',')) {
            y_read[row] = std::stod(token);
        }
        row++;
    }
    infile.close();

    cout << "hey" << endl;
    auto X_transpoze = Transpose(X_read, numRows, 14);
    cout << "oy" << endl;
    auto Xt_x = multiplyMatrices(X_transpoze, X_read, 14, numRows, 14);
    cout << "done" << endl;

    std::cout << "Features (X):" << std::endl;
    for (int i = 0; i < 14; ++i) {
        for (int j = 0; j < numCols - 1; ++j) {
            std::cout << Xt_x[i][j] << " ";
        }
        std::cout << std::endl;
    }

    delete[] X_read;
    delete[] y_read;

    return 0;
}


double* generateLaplaceNoise(double scale, int size) {
    std::random_device rd; //bu seed oluşturmak için
    std::mt19937 gen(rd()); //bu o seed'i kullanarak random bir sayı olışturuyor
    std::uniform_real_distribution<> dis(0.0, 1.0);//!!!!!!bu sıkıntı aslında  çünkü uint64 hali ile yapıyoruz işlemi

    double *noises = new double[size];
    for(int i=0; i<size; ++i){
        double u = dis(gen) - 0.5;//bu da o random sayıyı kullanıyor
        double noise = scale * ((u >= 0) ? -std::log(1 - 2*u) : std::log(1 + 2*u));
        noises[i] = noise;
    }
    return noises;
}



uint64_t *add_noise(Party *const proxy, uint64_t *numbers, uint32_t sz, uint32_t scale_32bit){
    const int length = int(sz);
    double scale = ConvertToDouble(uint64_t(scale_32bit));
    cout << "length is " << length << endl;
    cout << "scale is " << scale << endl;

    if (proxy->GetPRole() == helper){ //helper gürültünün imalinden mesul
        
        //bu kısımda gürültü imal edilmeli ve onun gizli payları oluşturulmalı.
        double *noise = new double[length];

        noise = generateLaplaceNoise(scale, length);
        
        uint64_t *noise_shares = ConvertToUint64(noise, length);

        uint64_t proxy1_share;
        uint64_t proxy2_share;
        
        unsigned char *ptr1 = proxy->GetBuffer1();
        unsigned char *ptr2 = proxy->GetBuffer2(); //HELPER'IN buffer'lara pointer'lar
        cout << "pointerlar falan tamam" << endl;
        for(int i=0; i<length; i++){ //bu optimal olmayabilir.
            cout << "for'a girdik" << endl;
            uint64_t tempShare = proxy->GenerateRandom();  //elimizde GenerateRandom da var, ikisi de aynı çalışıyor gibi.
            cout << "rastgele üretim" << endl;
            proxy1_share = tempShare; //ilki rastgele bir sayı olacak.
            proxy2_share = noise_shares[i] - tempShare; //ikincisi de o rastgele ile bizim göndermek istediğimizin farkı.

            AddValueToCharArray(proxy1_share, &ptr1); //Helper proxy1 için olan verisini kendi 1. buffer'ına ve
            AddValueToCharArray(proxy2_share, &ptr2); //proxy2'nin verisini kendi 2. buffer'ına yüklüyor
        }        
        cout << "for looptan çıktık" << endl;
        thread thr1 = thread(Send, proxy->GetSocketP1(), proxy->GetBuffer1(), length * 8);
        thread thr2 = thread(Send, proxy->GetSocketP2(), proxy->GetBuffer2(), length * 8);        

        thr1.join();
        thr2.join();

        delete[] noise;
        return nullptr;
    }
    else{
        cout << "add_noise proxy başla" << endl;
        uint64_t *noise_shares = new uint64_t[length];
        Receive(proxy->GetSocketHelper(), proxy->GetBuffer1(), 8*length); 
        cout << "add_noise proxy receive tamam" << endl;

        unsigned char *ptr = proxy->GetBuffer1(); //KENDİ birinci buffer'ı (birincinin buffer'ı değil!!)     

        for(int i=0; i<length; ++i){
            noise_shares[i] = ConvertToLong(&ptr);
        }
        cout << "add_noise proxy çevirmeler tamam" << endl;
        uint64_t *noiseful_shares = Add(proxy, numbers, noise_shares, length);


        delete[] noise_shares;
        return noiseful_shares;
    }
}










int main(){

}