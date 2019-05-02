// Utilities and system includes
#include <assert.h>
//#include <helper_string.h>  // helper for shared functions common to CUDA Samples
#include <iostream>
#include <fstream>
#include <string.h>
#include <iterator>
#include <thread>
#include <vector>
#include <omp.h>
#include <stdexcept>
#include <math.h>
// CUDA runtime
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cublasXt.h>
// CUDA and CUBLAS functions
#include <helper_functions.h>
#include <helper_cuda.h>




// Class variables
	static std::string outFile;
	static std::string distFile;
	static int threadCount = std::thread::hardware_concurrency() / 2;
	static double inflationParam = 3.0;
	static double epsilonParam = 0.001;
	static double zeroCutoff = 0.0001;
	static int strassenSize = 4;//1024;
	static bool isSim = false;
	static double** simMatrix;
	static double** simMatrixOriginal;
	static int nodeCount;
	std::vector<std::string> nodeList;

void normalizeMatrix(double** &mat) {
	double* nodeDegrees = new double[nodeCount];
	#pragma omp parallel for
	for(int i = 0; i < nodeCount; i++){
		nodeDegrees[i] = 0;
		for(int j = 0; j < nodeCount; j++) {
			nodeDegrees[i] += mat[i][j];
		}
	}
	#pragma omp parallel for
	for(int i = 0; i < nodeCount; i++){
		for(int j = 0; j < nodeCount; j++) {
			mat[i][j] /= nodeDegrees[i];
		}
	}
}
void setEntriesZero(double** &mat) {
	#pragma omp parallel for
	for(int i = 0; i < nodeCount; i++){
		for(int j = 0; j < nodeCount; j++) {
			if(mat[i][j] < zeroCutoff) {
				mat[i][j] = 0;
			}
		}
	}
}
double calcMatNorm(double** matA, double** matB) {
	double sum = 0;
	for(int i = 0; i < nodeCount; i++) {
		for(int j = 0; j < nodeCount; j++) {
			sum += pow(matA[i][j]-matB[i][j], 2);
		}
	}
	return sqrt(sum);
}
void calcMatInflation(double** &mat) {
	#pragma omp parallel for
	for(int i = 0; i < nodeCount; i++){
		for(int j = 0; j < nodeCount; j++) {
			mat[i][j] = pow(mat[i][j], inflationParam);
		}
	}
	normalizeMatrix(mat);
}
void matToArray(double** &mat, int matrixSize, double* &array){
	array = new double[matrixSize*matrixSize];
	#pragma omp parallel for
	for(int i = 0; i < matrixSize; i++){
		for(int j = 0; j < matrixSize; j++){
			array[i*matrixSize+j] = mat[i][j]; 
		}	
	}
}
void arrayToMat(double* array, int matrixSize, double** &mat){
	mat = new double*[matrixSize];
	#pragma omp parallel for
	for(int i = 0; i < matrixSize; i++){
		mat[i] = new double[matrixSize];
		for(int j = 0; j < matrixSize; j++){
			mat[i][j] = array[i*matrixSize+j]; 
		}	
	}
}
std::vector< std::pair< std::string,std::vector< std::string > > > findClusters(double** mat){
	std::vector< std::pair< std::string,std::vector< std::string > > > clusterList;
	std::vector<int> attractorList;
	for(int i = 0; i < nodeCount; i++){
		if(mat[i][i] > 0){
			attractorList.push_back(i);
			clusterList.push_back(std::make_pair(nodeList[i], std::vector< std::string > ()));
		}
	}
	#pragma omp parallel for
	for(int i = 0; i < attractorList.size(); i++){
		for(int j = 0; j < nodeCount; j++) {
			if(mat[j][attractorList[i]] > 0) {
					clusterList[i].second.push_back(nodeList[j]);
			}
		}
	}
	return clusterList;	
}
int findPos(std::vector<std::string> vec, std::string element){
	return std::distance(vec.begin(), find(vec.begin(), vec.end(), element));
}
void analyzeCluster(std::pair<std::string,std::vector<std::string>> cluster, std::ofstream &out){
	std::cout << "Cluster-Size = " << cluster.second.size() << std::endl;
	int posFirst = findPos(nodeList, cluster.first);
	std::vector<int> posSecond;
	for(int i = 0; i < cluster.second.size(); i++){
		posSecond.push_back(findPos(nodeList, cluster.second[i]));
	}
	out << "Cluster " << cluster.first << " consists of " << cluster.second.size() << " entries of " << nodeCount << " elements in total\n";
	double meanPairwiseSimOriginal = 0;
	double meanRepSimOriginal = 0;
	double meanPairwiseSimCluster = 0;
	double meanRepSimCluster = 0;	
	#pragma omp parallel for reduction(+:meanRepSimOriginal,meanRepSimCluster)
	for(int i = 0; i < cluster.second.size(); i++){
		meanRepSimOriginal += simMatrixOriginal[posFirst][posSecond[i]];
		meanRepSimCluster += simMatrix[posFirst][posSecond[i]];
	}
	#pragma omp parallel for reduction(+:meanPairwiseSimOriginal,meanPairwiseSimCluster)
	for(int i = 0; i < cluster.second.size(); i++){
		double localSumOriginal = 0.0;
		double localSumCluster = 0.0;
		for(int j = 0; j < cluster.second.size(); j++){
			localSumOriginal += simMatrixOriginal[posSecond[i]][posSecond[j]];
			localSumCluster += simMatrix[posSecond[i]][posSecond[j]];
		}
		meanPairwiseSimOriginal += localSumOriginal;
		meanPairwiseSimCluster += localSumCluster;	
	}
	meanRepSimOriginal /= cluster.second.size();
	meanRepSimCluster /= cluster.second.size();
	meanPairwiseSimOriginal /= (cluster.second.size()*cluster.second.size());
	meanPairwiseSimCluster /= (cluster.second.size()*cluster.second.size());
	out << "Mean Pairwise-Sim before clustering = " << meanPairwiseSimOriginal << "\n";
	out << "Mean Pairwise-Sim after clustering = " << meanPairwiseSimCluster << "\n";
	out << "Mean Sim-to-Rep before clustering = " << meanRepSimOriginal << "\n";
	out << "Mean Sim-to-Rep after clustering = " << meanRepSimCluster << "\n";
	out << "Rep = " << cluster.first << "\n";
	out << "[" << cluster.second[0];
	for(int i = 1; i < cluster.second.size(); i++){
		out << ", " << cluster.second[i];
	}
	out << "]\n";
}
////////////////////////////////////////////////////////////////////////////////
//! Run matrix multiplication using CUBLASXt
////////////////////////////////////////////////////////////////////////////////
int matrixMult_Cuda_CublasXt(double** matA, double** matB, double** &matC, int matrixSize){
	double* arrayA;
	double* arrayB;
	double* arrayC;
	cublasStatus_t status;
	//Convert matrices to arrays
	matToArray(matA, matrixSize, arrayA);
	matToArray(matB, matrixSize, arrayB);
	arrayC = new double[matrixSize*matrixSize];
	//Create a CULBASXt handle
	cublasXtHandle_t handle;
   	status = cublasXtCreate(&handle);
	if (status != CUBLAS_STATUS_SUCCESS) {
    	printf("!!!! CUBLASXT initialization error\n");
    	return EXIT_FAILURE;
	}
	//Select GPUs
	int *devices = NULL;
	int num_of_devices = 0;
	checkCudaErrors(cudaGetDeviceCount(&num_of_devices));
	devices = new int[num_of_devices];
	for(int i = 0; i < num_of_devices; i++){
		devices[i] = i;
	}
	status = cublasXtDeviceSelect(handle, num_of_devices, devices);
	if (status != CUBLAS_STATUS_SUCCESS) {
    	printf("!!!! CUBLASXT device selection error\n");
    	return EXIT_FAILURE;
	}
	// Execute dgemm
	const double alpha = 1.0;
	const double beta = 0.0;
	status = cublasXtDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, matrixSize, matrixSize, matrixSize, &alpha, arrayA, matrixSize, arrayB, matrixSize, &beta, arrayC, matrixSize);
	if (status != CUBLAS_STATUS_SUCCESS) {
    	printf("!!!! CUBLASXT kernel execution error\n");
    	return EXIT_FAILURE;
	}
	//Clean up
	status = cublasXtDestroy(handle);
	if (status != CUBLAS_STATUS_SUCCESS) {
    	printf("!!!! CUBLASXT shutdown error\n");
    	return EXIT_FAILURE;
	}
	//Convert array to matrix
	arrayToMat(arrayC, matrixSize, matC);
	return EXIT_SUCCESS;
}

void readMatrix(){
	std::string line;
	std::ifstream in;
 	in.open(distFile);
	getline(in,line);
	std::istringstream iss(line);
	std::vector<std::string> tmpList(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());
	nodeCount = tmpList.size()-1;
	std::cout <<"Reading matrix with " << nodeCount << " nodes\n";
	simMatrix = new double*[nodeCount];
	for(int i = 0; i < nodeCount; i++){
		simMatrix[i] = new double[nodeCount];
	}
	double* nodeDegrees = new double[nodeCount];
	for(int i = 0; i < nodeCount; i++){
		for(int j = 0; j < nodeCount; j++){
			simMatrix[i][j] = 0.0;
		}
		nodeDegrees[i] = 0.0;
	}
	std::cout <<"Matrix init start\n";
	double maxValue = 0; 
	for(int i = 0; i < nodeCount; i++){
		nodeList.push_back(tmpList[0]);
		for(int j = i+1; j < nodeCount; j++){
			double tmp = std::stod(tmpList[j+1]);
			simMatrix[i][j] = tmp;
			simMatrix[j][i] = tmp;
			nodeDegrees[i] += tmp;
			nodeDegrees[j] += tmp;
			if(tmp > maxValue) {
				maxValue = tmp;
			}
		}
		getline(in,line);
		iss = std::istringstream(line);
	 	tmpList = std::vector<std::string>(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());
	}
	in.close();

	//Convert distances into similarities
	if(!isSim) { 
			
	}
	//Add loops to the matrix
	for(int i = 0; i < nodeCount; i++) {
		simMatrix[i][i] = 1;
		nodeDegrees[i] += 1;
	}

	//Remove entries smaller than cutoff
	#pragma omp parallel for
	for(int i = 0; i < nodeCount; i++){
		for(int j = 0; j < nodeCount; j++) {
			if(simMatrix[i][j] < 0.2) {
				simMatrix[i][j] = 0;
			}
		}
	}

	//Copy original matrix
	simMatrixOriginal = new double*[nodeCount];
	#pragma omp parallel for
	for(int i = 0; i < nodeCount; i++){
		simMatrixOriginal[i] = new double[nodeCount];
		for(int j = 0; j < nodeCount; j++) {
			simMatrixOriginal[i][j] = simMatrix[i][j];
		}
	}

	//Normalizing matrix
	normalizeMatrix(simMatrix);
	std::cout << "Matrix init done\n"; 
}
void parseArgs(int argc, char **argv){
	for(int i = 1; i < argc-1; i++){
		std::string arg = argv[i];
		if(arg == "-out"){
			outFile = argv[i+1];
			i++;
		}
		else if(arg == "-threads"){
			threadCount = std::stoi(argv[i+1]);
			i++;
		}
		else if(arg == "-sim"){
			isSim = true;
		}
		else if(arg == "-epsilon"){
			epsilonParam = std::stod(argv[i+1]);
			i++;
		}
		else if(arg == "-inflation"){
			inflationParam = std::stod(argv[i+1]);
			i++;
		}
		else if(arg == "-zero"){
			zeroCutoff = std::stod(argv[i+1]);
			i++;
		}
		else{
			std::cout << "Unknown argument " << argv[i] << "\n";
		}
	}
	distFile = argv[argc-1];
	if(outFile.empty()){
		outFile = distFile + ".cls";	
	}
	if(epsilonParam < 0) {
		throw std::invalid_argument("Epsilon needs to be greater than 0");
	}
	if(inflationParam < 1) {
		throw std::invalid_argument("Inflation parameter needs to be at least 1");
	}
	if(zeroCutoff < 0) {
		throw std::invalid_argument("Cutoff value needs to be greater than 0");
	}
}
void printHelp(){
	std::cout << "Usage: ./Markov_Clustering_Cuda {Options} Dist.Matrix\n";
	std::cout << "Options:\n";
	std::cout << "-out file\t name of outputfile (default Dist.Matrix.cls)\n";
	std::cout << "-threads number\t number of threads to use (default = Number_of_Cores / 2)\n";
	std::cout << "-sim\t indicates that the matrix contains similarity values\n";
	std::cout << "-epsilon number\t number to use as stop criterium (default = 0.001)\n";
	std::cout << "-inflation number\t number to use as inflation parameter\n";
	std::cout << "-zero number\t cutoff-value for removing connections (default = 0.0001)\n";
}
////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv){
	if(argc < 2){
		printHelp();
		return 0;
	}
	try{
		parseArgs(argc, argv);
	}
	catch(const std::invalid_argument& e){
		std::cout << "Could not parse arguments\n";
		std::cout << e.what() << "\n";
		return 0;
	}
	omp_set_dynamic(0);
	omp_set_num_threads(threadCount);
	readMatrix();
	double** oldMatrix = new double*[nodeCount];
	for(int i = 0; i < nodeCount; i++){
		oldMatrix[i] = new double[nodeCount];
	}
	int iteration = 0;
	do{
		std::cout << "Starting iteration " << iteration << "\n";
		#pragma omp parallel for
		for(int i = 0; i < nodeCount; i++){
			for(int j = 0; j < nodeCount; j++){
				oldMatrix[i][j] = simMatrix[i][j];
			}
		}
		std::cout << "\tExpand\n";
		int result = matrixMult_Cuda_CublasXt(simMatrix, simMatrix, simMatrix, nodeCount);
		if(result != EXIT_SUCCESS){
			std::cout << "An error occured during the matrix multiplication\n";
			return -1;
		}
		std::cout << "\tInflate\n";
		calcMatInflation(simMatrix);
		setEntriesZero(simMatrix);
		std::cout <<"Finished iteration " << iteration << "\n";
		iteration++;
	} while(calcMatNorm(simMatrix,oldMatrix) > epsilonParam);
	for(int i = 0; i < nodeCount; i++) {
		double sum = 0; 
		for(int j = 0; j < nodeCount; j++) {
			sum += simMatrix[i][j];
		}
		if(abs(sum - 1) > 0.00001) {
			std::cout << "Difference in row-sum for entry " << i << " = " << nodeList[i] << " Sum = " << sum << "\n";
		}
	}
	std::cout << "Finding Clusters\n";
	std::vector< std::pair< std::string,std::vector< std::string > > > clusterList = findClusters(simMatrix);
	std::cout << "Printing Results\n";
	std::ofstream out(outFile);
	out << "Cluster-Ergebnisse:\n";
	for(int i =0; i < clusterList.size(); i++){
		std::cout << "Cluster " << i << " of " << clusterList.size() << "\n";
		analyzeCluster(clusterList[i], out);
	}
	out.close();
	std::cout << "Finished\n";
	return 0;
}

