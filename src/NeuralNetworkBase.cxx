#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "UHH2/BoostedHiggsToWW/include/NeuralNetworkBase.hpp"

using namespace std;

NeuralNetworkBase::NeuralNetworkBase(const string& path_) : path(path_) {
  LoadLayers();
  extractVariables();
  extractNormalization();
  // PrintVariables();
  // PrintNormalization();
};


void NeuralNetworkBase::extractVariables() {
  string line;
  ifstream myfile(finfo, ios::in);
  while (!myfile.eof()) {
    std::getline(myfile, line);
    if (line.find("variables") != std::string::npos) {
      line = line.substr(line.find("[")+1, line.find("]")-line.find("[")-1);
      TString tok; int from = 0;
      while (((TString)line).Tokenize(tok, from, ", ")) {
        variables.push_back(tok(1,tok.Length()-2));
      }
    }
  }
  myfile.close();
}



void NeuralNetworkBase::extractNormalization() {
  string line;
  ifstream myfile(fnorm, ios::in);
  std::getline(myfile, line);
  while (std::getline(myfile, line) && line!="") {
    vector<TString> temp;
    TString tok; int from = 0;
    while (((TString)line).Tokenize(tok, from, " ")) {
      temp.push_back(tok);
    }
    Normalizer.insert(std::make_pair(temp[0], std::make_pair(temp[1],std::make_pair(temp[2].Atof(), temp[3].Atof()))));
    TString input = temp[0];
  }
  Normalizer.insert(std::make_pair("jetBtag", std::make_pair("MinMaxScaler",std::make_pair(0., 1.))));
  myfile.close();


}


void NeuralNetworkBase::PrintVariables(){

  for (size_t i = 0; i < variables.size(); i++) {
    std::cout << variables[i] << '\t';
  }
  std::cout << std::endl;

}

void NeuralNetworkBase::PrintNormalization() {
  std::cout << "Normalizer" << '\n';
  for(const auto & input : variables){
    std::cout << input << "\t" << Normalizer[input].first << "\t" << Normalizer[input].second.first << "\t" << Normalizer[input].second.second << '\n';
  }
  std::cout << std::endl;
}



void NeuralNetworkBase::extractLayers() {
  string line;
  vector<TString> temp;
  ifstream myfile(finfo, ios::in);
  while (!myfile.eof()) {
    std::getline(myfile, line);
    if (line.find("layers") != std::string::npos) {
      line = line.substr(line.find("[")+1, line.find("]")-line.find("[")-1);
      TString tok; int from = 0;
      while (((TString)line).Tokenize(tok, from, ", ")) {
        layers_dim.push_back(tok.Atoi());
      }
    }
    if (line.find("sample_names") != std::string::npos) {
      line = line.substr(line.find("[")+1, line.find("]")-line.find("[")-1);
      TString tok; int from = 0;
      while (((TString)line).Tokenize(tok, from, ", ")) {
        temp.push_back(tok);
      }
    }
  }
  layers_dim.push_back(temp.size());
  myfile.close();
}


vector< float > NeuralNetworkBase::readLineWeights(TString token) {
  vector< float > temp;
  TString tok; int from = 0;
  while (token.Tokenize(tok, from, ",\t")) {
    temp.push_back(tok.Atof());
  }
  return temp;
}







Matrix2D NeuralNetworkBase::MatrixMultiplication( Matrix2D A, Matrix2D B){
  if (A[0].size()!=B.size()) throw std::invalid_argument( "received invalid argument" );
  Matrix2D C (A.size(), vector<float> (B[0].size(),0));
  for (size_t i = 0; i < C.size(); i++) {
    for (size_t j = 0; j < C[0].size(); j++) {
      for (size_t k = 0; k < B.size(); k++) {
        C[i][j] += A[i][k]*B[k][j];
      }
    }
  }
  return C;
}

Matrix2D NeuralNetworkBase::MatrixAddition( Matrix2D A, Matrix2D B){
  if (A.size()!=B.size() || A[0].size()!=B[0].size() ) throw std::invalid_argument( "received invalid argument" );
  Matrix2D C (A.size(), vector<float> (A[0].size(),0));
  for (size_t i = 0; i < C.size(); i++) {
    for (size_t j = 0; j < C[0].size(); j++) {
      C[i][j] = A[i][j]+B[i][j];
    }
  }
  return C;
}

Matrix2D NeuralNetworkBase::ApplyActivation( Matrix2D A, TString activation){
  if (activation=="ID") {
    return A;
  }
  Matrix2D C (A.size(), vector<float> (A[0].size(),0));
  if (activation=="relu") {
    for (size_t i = 0; i < C.size(); i++) {
      for (size_t j = 0; j < C[0].size(); j++) {
        C[i][j] = (A[i][j]>0)?A[i][j]:0;
      }
    }
  }
  if (activation=="softmax") {
    for (size_t i = 0; i < C.size(); i++) {
      double sum = 0.;
      for (size_t j = 0; j < C[0].size(); j++) sum += exp(A[i][j]);
      for (size_t j = 0; j < C[0].size(); j++) {
        C[i][j] = exp(A[i][j])/sum;
      }
    }
  }
  return C;
}



void NeuralNetworkBase::LoadLayers() {
  string line;
  vector< string > lines;
  ifstream myfile(fname, ios::in);
  while (std::getline(myfile, line)) lines.push_back(line);
  myfile.close();

  for (size_t i = 0; i < lines.size(); i++) {
    line = lines[i];

    if (line.find("New Layer") != std::string::npos) {
      line = lines[i+=1];
      layers_name.push_back((TString)line);
    }

    if (line.find("weights") != std::string::npos) {
      vector< vector< float > > weights;
      while (true) {
        line = lines[i+=1];
        if (line.find("bias") != std::string::npos) {
          break;
        }
        line = line.substr(0,line.size()-1);
        // vector< float > temp;
        // readLineWeights((TString)line,temp);
        weights.push_back(readLineWeights((TString)line));
      }
      layers.push_back(weights);
    }

    if (line.find("bias") != std::string::npos) {
      line = lines[i+=1];
      line = line.substr(0,line.size()-1);
      vector< vector< float > > bias;
      bias.push_back(readLineWeights((TString)line));
      biases.push_back(bias);
    }

    if (line.find("activation") != std::string::npos) {
      line = line.substr(line.find("\t"),line.size());
      line.erase(std::remove(line.begin(), line.end(), '\t'), line.end());
      activations.push_back((TString)line);
    }

    if (line.find("epsilon") != std::string::npos) {
      vector< float > epsilon;
      line = lines[i+=1];
      epsilon.push_back(((TString)line).Atof());
      line = lines[i+=2];
      vector< float > gamma = readLineWeights((TString)line);
      line = lines[i+=2];
      vector< float > beta = readLineWeights((TString)line);
      line = lines[i+=2];
      vector< float > moving_mean = readLineWeights((TString)line);
      line = lines[i+=2];
      vector< float > moving_variance = readLineWeights((TString)line);
      vector< vector< float > > matrix(moving_variance.size(), vector< float > (moving_variance.size(), 0.) );
      vector< vector< float > > bias(1, vector< float > (moving_variance.size(), 0.) );
      for (size_t i = 0; i < moving_variance.size(); i++) {
        matrix[i][i] = gamma[i]/sqrt(moving_variance[i]+epsilon[0]);
        bias[0][i] = - moving_mean[i]*gamma[i]/sqrt(moving_variance[i]+epsilon[0])+beta[i];
      }
      layers.push_back(matrix);
      biases.push_back(bias);
      activations.push_back((TString)("ID"));
    }
  }
}

void NeuralNetworkBase::CreateInputs(const TopJet& jet) {
  std::vector<float> v;
  for(const auto & var : variables){
    double val;
    if (var == "jetPt") val = jet.pt();
    else if (var == "jetEta") val = jet.eta();
    else if (var == "jetPhi") val = jet.phi();
    else if (var == "jetMass") val = jet.v4().M();
    else if (var == "jetEnergy") val = jet.energy();
    else if (var == "jetTau1") val = jet.tau1();
    else if (var == "jetTau2") val = jet.tau2();
    else if (var == "jetTau3") val = jet.tau3();
    else if (var == "jetTau4") val = jet.tau4();
    else if (var == "ncandidates") val = jet.chargedMultiplicity()+jet.neutralMultiplicity();
    else if (var == "jetBtag") val = jet.btag_combinedSecondaryVertex();
    else if (var == "jetTau21") val = (jet.tau1()!=0) ? jet.tau2()/jet.tau1() : 0;
    else if (var == "jetTau31") val = (jet.tau1()!=0) ? jet.tau3()/jet.tau1() : 0;
    else if (var == "jetTau41") val = (jet.tau1()!=0) ? jet.tau4()/jet.tau1() : 0;
    else if (var == "jetTau32") val = (jet.tau2()!=0) ? jet.tau3()/jet.tau2() : 0;
    else if (var == "jetTau42") val = (jet.tau2()!=0) ? jet.tau4()/jet.tau2() : 0;
    else if (var == "jetTau43") val = (jet.tau3()!=0) ? jet.tau4()/jet.tau3() : 0;
    else throw std::invalid_argument( "received invalid argument" );
    v.push_back(val);
  }
  inputs.clear();
  inputs.push_back(v);
}

void NeuralNetworkBase::NormalizeInput() {
  for (size_t i = 0; i < variables.size(); i++) {
    if (Normalizer[variables[i]].first == "StandardScaler" || Normalizer[variables[i]].first == "StandardScalerNoMean") {
      inputs[0][i] = (inputs[0][i]+ Normalizer[variables[i]].second.first) * Normalizer[variables[i]].second.second;
    } else if (Normalizer[variables[i]].first == "MinMaxScaler") {
      inputs[0][i] = inputs[0][i]* Normalizer[variables[i]].second.second + Normalizer[variables[i]].second.first;
    } else {
      throw std::invalid_argument( Normalizer[variables[i]].first+" Scaler is not Implemented");
    }
  }
}

void NeuralNetworkBase::Apply(const TopJet& jet) {
  CreateInputs(jet);
  // std::cout << "Input" << '\n';
  // for (size_t i = 0; i < inputs[0].size(); i++) {
  //   std::cout << inputs[0][i] << '\t';
  // }
  // std::cout << '\n';
  NormalizeInput();
  outputs.clear();
  // outputs = inputs;
  outputs = inputs;
  // std::cout << "NORM" << '\n';
  // for (size_t i = 0; i < outputs[0].size(); i++) {
  //   std::cout << outputs[0][i] << '\t';
  // }
  // std::cout << '\n';
  for (size_t i = 0; i < layers_name.size(); i++) {
    outputs = MatrixMultiplication(outputs,layers[i]);
    outputs = MatrixAddition(outputs,biases[i]);
    outputs = ApplyActivation(outputs,activations[i]);
  }
  std::cout << "Output" << '\n';
  for (size_t i = 0; i < outputs[0].size(); i++) {
    std::cout << outputs[0][i] << '\t';
  }
  std::cout << '\n';
}
