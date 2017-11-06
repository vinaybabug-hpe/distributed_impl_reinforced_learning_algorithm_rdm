// Copyright 2009 Google Inc.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <ctime>
#include <iostream>
#include "pspectral/common.h"
#include "pspectral/compute_distance.h"


namespace learning_psc {
Document::Document() : two_norm_sq(-1) {
}
void Document::Encode(string* s) {
  s->clear();
  IntToString(iv.size(), s);
  for (int i = 0; i < iv.size(); ++i) {
    IntToString(iv[i].index, s);
    DoubleToString(iv[i].value, s);
  }
  DoubleToString(two_norm_sq, s);
}
void Document::Decode(const string& s) {
  iv.clear();
  const char* pointer = s.c_str();
  int iv_size = StringToInt(pointer, sizeof(int));
  iv.resize(iv_size);
  pointer += sizeof(int);
  for (int i = 0; i < iv_size; ++i) {
    iv[i].index = StringToInt(pointer, sizeof(int));
    pointer += sizeof(int);
    iv[i].value = StringToDouble(pointer, sizeof(double));
    pointer += sizeof(double);
  }
  two_norm_sq = StringToDouble(pointer, sizeof(double));
  pointer += sizeof(double);
  CHECK_EQ(*pointer, '\0');
}

void Document::ComputeTwoNorm() {
  two_norm_sq = 0;
  for (int i = 0; i < iv.size(); ++i) {
    two_norm_sq += iv[i].value * iv[i].value;
  }
}

ComputeDistance::ComputeDistance(int t_nearest_neighbor, MPI_Comm new_comm)
  : t_nearest_neighbor_(t_nearest_neighbor),
	communicator (new_comm) {
  MPI_Comm_rank(communicator, &myid_);
  MPI_Comm_size(communicator, &pnum_);
  d_printf("myid_=%d pnum_=%d\n", myid_, pnum_);
}
/*void ComputeDistance::Read(const string& filename) {
  ifstream fin2(filename.c_str());
  string line;
  int i = 0;
  while (getline(fin2, line)) {  // Each line is a training document.
    if (line.size() > 0 &&      // Skip empty lines.
        line[0] != '\r' &&      // Skip empty lines.
        line[0] != '\n' &&      // Skip empty lines.
        line[0] != '#') {       // Skip comment lines.


    	d_printf(" READ: %s\n", line.c_str());

      if (i++ % pnum_ != myid_) {
        continue;
      }


      Document document;
      istringstream ss(line);
      int index;
      double value;
      char colon;
      while (ss >> index >> colon >> value) {  // Load and init a document.
        document.iv.push_back(IndexValue(index, value));
      }
      document.ComputeTwoNorm();
      docs_.push_back(document);
    }
  }
  num_total_docs_ = i;

  fin2.close();

}*/

void ComputeDistance::Read(char* dataInputBuffer) {
  std::istringstream fin2(dataInputBuffer);
  string line;
  int i = 0;
  while (getline(fin2, line)) {  // Each line is a training document.
    if (line.size() > 0 &&      // Skip empty lines.
        line[0] != '\r' &&      // Skip empty lines.
        line[0] != '\n' &&      // Skip empty lines.
        line[0] != '#') {       // Skip comment lines.


//    	d_printf(" READ: %s\n", line.c_str());

      if (i++ % pnum_ != myid_) {
        continue;
      }


      Document document;
      istringstream ss(line);
      int index;
      double value;
      char colon;
      while (ss >> index >> colon >> value) {  // Load and init a document.
        document.iv.push_back(IndexValue(index, value));
      }
      document.ComputeTwoNorm();
      docs_.push_back(document);
    }
  }
  num_total_docs_ = i;


}

void ComputeDistance::ParallelComputeEuclidean() {


  vector<TopN<IndexValue, CmpIndexValueByValueAsc>*> topn;
  for (int j = 0; j < docs_.size(); ++j) {
    topn.push_back(new TopN<IndexValue, CmpIndexValueByValueAsc>(t_nearest_neighbor_));
  }


  time_t t = time(NULL);
  for (int i = 0; i < num_total_docs_; ++i) {


    // Where this docs is stored.
    int which_computer = i % pnum_;
    // Broadcast this document to all computer
    Document doc;
    if (which_computer == myid_) {
      doc = docs_[i / pnum_];
    }



    BroadcastDocument(&doc, which_computer);


    for (int j = 0; j < docs_.size(); ++j) {
      if (which_computer == myid_ && i / pnum_ == j) {
        // Do not compute myself to myself
        continue;
      }
      double distance = sqrt(doc.two_norm_sq +
             docs_[j].two_norm_sq -
             2 * InnerProduct(doc, docs_[j]));
      topn[j]->Insert(IndexValue(i, distance));

//      d_printf("[@%d] DISTANCE = %f\n", myid_, distance);

    }



    if (myid_ == 0 && i % 1000 == 0 && PRINT_DEBUG) {
      time_t now = time(NULL);
      LOG(INFO) << "Progress: " << i << " / " << num_total_docs_ << " time_elapsed:"
                << now - t << "s" << std::endl;
    }


  }
  distance_.clear();
  distance_.resize(topn.size());
  for (int i = 0; i < topn.size(); ++i) {
    topn[i]->Extract(&distance_[i]);
    sort(distance_[i].begin(), distance_[i].end(), CmpIndexValueByIndexAsc());
  }
  for (int j = 0; j < docs_.size(); ++j) {
    delete topn[j];
  }
}

void ComputeDistance::ParallelComputeManhattan() {
  vector<TopN<IndexValue, CmpIndexValueByValueAsc>*> topn;
  for (int j = 0; j < docs_.size(); ++j) {
    topn.push_back(new TopN<IndexValue, CmpIndexValueByValueAsc>(t_nearest_neighbor_));
  }
  time_t t = time(NULL);
  for (int i = 0; i < num_total_docs_; ++i) {
    // Where this docs is stored.
    int which_computer = i % pnum_;
    // Broadcast this document to all computer
    Document doc;
    if (which_computer == myid_) {
      doc = docs_[i / pnum_];
    }
    BroadcastDocument(&doc, which_computer);
    for (int j = 0; j < docs_.size(); ++j) {
      if (which_computer == myid_ && i / pnum_ == j) {
        // Do not compute myself to myself
        continue;
      }
      double distance = cityDistance(doc, docs_[j]);
      topn[j]->Insert(IndexValue(i, distance));
    }
    if (myid_ == 0 && i % 1000 == 0 && PRINT_DEBUG) {
      time_t now = time(NULL);
      LOG(INFO) << "Progress: " << i << " / " << num_total_docs_ << " time_elapsed:"
                << now - t << "s" << std::endl;
    }
  }
  distance_.clear();
  distance_.resize(topn.size());
  for (int i = 0; i < topn.size(); ++i) {
    topn[i]->Extract(&distance_[i]);
    sort(distance_[i].begin(), distance_[i].end(), CmpIndexValueByIndexAsc());
  }
  for (int j = 0; j < docs_.size(); ++j) {
    delete topn[j];
  }
}

void ComputeDistance::ParallelComputeCorrelation() {
  vector<TopN<IndexValue, CmpIndexValueByValueAsc>*> topn;
  for (int j = 0; j < docs_.size(); ++j) {
    topn.push_back(new TopN<IndexValue, CmpIndexValueByValueAsc>(t_nearest_neighbor_));
  }
  time_t t = time(NULL);
  for (int i = 0; i < num_total_docs_; ++i) {
    // Where this docs is stored.
    int which_computer = i % pnum_;
    // Broadcast this document to all computer
    Document doc;
    if (which_computer == myid_) {
      doc = docs_[i / pnum_];
    }
    BroadcastDocument(&doc, which_computer);
    for (int j = 0; j < docs_.size(); ++j) {
      if (which_computer == myid_ && i / pnum_ == j) {
        // Do not compute myself to myself
        continue;
      }
      double sum1 = calcSummation(doc);
      double sum2 = calcSummation(docs_[j]);
      double result = InnerProduct(doc, docs_[j]);
      double denom1 = doc.two_norm_sq;
      double denom2 = docs_[j].two_norm_sq;

      result -= sum1 * sum2;
      denom1 -= sum1 * sum1;
      denom2 -= sum2 * sum2;

      result = result / sqrt(denom1*denom2);
      result = 1. - result;

      double distance = result;
      topn[j]->Insert(IndexValue(i, distance));
    }
    if (myid_ == 0 && i % 1000 == 0 && PRINT_DEBUG) {
      time_t now = time(NULL);
      LOG(INFO) << "Progress: " << i << " / " << num_total_docs_ << " time_elapsed:"
                << now - t << "s" << std::endl;
    }
  }
  distance_.clear();
  distance_.resize(topn.size());
  for (int i = 0; i < topn.size(); ++i) {
    topn[i]->Extract(&distance_[i]);
    sort(distance_[i].begin(), distance_[i].end(), CmpIndexValueByIndexAsc());
  }
  for (int j = 0; j < docs_.size(); ++j) {
    delete topn[j];
  }
}

void ComputeDistance::ParallelComputeCosine() {
  vector<TopN<IndexValue, CmpIndexValueByValueAsc>*> topn;
  for (int j = 0; j < docs_.size(); ++j) {
    topn.push_back(new TopN<IndexValue, CmpIndexValueByValueAsc>(t_nearest_neighbor_));
  }
  time_t t = time(NULL);
  for (int i = 0; i < num_total_docs_; ++i) {
    // Where this docs is stored.
    int which_computer = i % pnum_;
    // Broadcast this document to all computer
    Document doc;
    if (which_computer == myid_) {
      doc = docs_[i / pnum_];
    }
    BroadcastDocument(&doc, which_computer);
    for (int j = 0; j < docs_.size(); ++j) {
      if (which_computer == myid_ && i / pnum_ == j) {
        // Do not compute myself to myself
        continue;
      }

      double denom1 = doc.two_norm_sq;
      double denom2 = docs_[j].two_norm_sq;
      double numerator = InnerProduct(doc, docs_[j]);
      double distance =  numerator / (sqrt(denom1)*sqrt(denom2));

      topn[j]->Insert(IndexValue(i, distance));
    }
    if (myid_ == 0 && i % 1000 == 0 && PRINT_DEBUG) {
      time_t now = time(NULL);
      LOG(INFO) << "Progress: " << i << " / " << num_total_docs_ << " time_elapsed:"
                << now - t << "s" << std::endl;
    }
  }
  distance_.clear();
  distance_.resize(topn.size());
  for (int i = 0; i < topn.size(); ++i) {
    topn[i]->Extract(&distance_[i]);
    sort(distance_[i].begin(), distance_[i].end(), CmpIndexValueByIndexAsc());
  }
  for (int j = 0; j < docs_.size(); ++j) {
    delete topn[j];
  }
}

/*
void ComputeDistance::Write(const string& file) {
  for (int i = 0; i < num_total_docs_; ++i) {
    if (i % pnum_ == myid_) {
      int local_index = i / pnum_;
      std::ofstream fout;
      if (i == 0) {
        fout.open(file.c_str());
      } else {
        fout.open(file.c_str(), std::ios::app);
      }
      for (int j = 0; j < distance_[local_index].size(); ++j) {
        if (j != 0) {
          fout << " ";
          std::cout << " ";
        }
        std::cout << distance_[local_index][j].index << ":" << distance_[local_index][j].value;
        fout << distance_[local_index][j].index << ":" << distance_[local_index][j].value;
      }
      std::cout << std::endl;
      fout << std::endl;
      fout.close();
    }
    MPI_Barrier(communicator);
  }
}
*/

char* ComputeDistance::Write() {

  std::stringstream dataBuffer;
  const char *dataPtr;

  for (int i = 0; i < num_total_docs_; ++i) {
    if (i % pnum_ == myid_) {
      int local_index = i / pnum_;
//      std::ofstream fout;
      if (i == 0) {
//        fout.open(file.c_str());
      } else {
//        fout.open(file.c_str(), std::ios::app);
      }
      for (int j = 0; j < distance_[local_index].size(); ++j) {
        if (j != 0) {
        	dataBuffer << " ";
//          std::cout << " ";
        }
//        std::cout << distance_[local_index][j].index << ":" << distance_[local_index][j].value;
        dataBuffer << distance_[local_index][j].index << ":" << distance_[local_index][j].value;
      }
//      std::cout << std::endl;
      dataBuffer << std::endl;
//      fout.close();
    }
    MPI_Barrier(communicator);
  }

  dataPtr = new char[dataBuffer.str().size()+1];
  strcpy(dataPtr, dataBuffer.str().c_str());
  return dataPtr;
}


void ComputeDistance::BroadcastDocument(Document* doc, int root) {
  string s;
  if (myid_ == root) {
    doc->Encode(&s);
  }
  int s_size = s.size();
  MPI_Bcast(&s_size, 1, MPI_INT, root, communicator);
  if (myid_ != root) {
    s.resize(s_size);
  }
  MPI_Bcast(&s[0], s_size, MPI_CHAR, root, communicator);
  doc->Decode(s);
}

double ComputeDistance::InnerProduct(const Document& doc1, const Document& doc2) {
  double inner_product = 0;
  int it1 = 0;
  int it2 = 0;
  while (it1 < doc1.iv.size() && it2 < doc2.iv.size()) {
    if (doc1.iv[it1].index == doc2.iv[it2].index) {
      inner_product += doc1.iv[it1].value * doc2.iv[it2].value;
      ++it1;
      ++it2;
    } else if (doc1.iv[it1].index < doc2.iv[it2].index) {
      ++it1;
    } else {
      ++it2;
    }
  }
  return inner_product;
}

double ComputeDistance::cityDistance(const Document& doc1, const Document& doc2) {
  double result = 0;
  int it1 = 0;
  int it2 = 0;
  while (it1 < doc1.iv.size() && it2 < doc2.iv.size()) {
    if (doc1.iv[it1].index == doc2.iv[it2].index) {
    	result += fabs(doc1.iv[it1].value - doc2.iv[it2].value);
      ++it1;
      ++it2;
    } else if (doc1.iv[it1].index < doc2.iv[it2].index) {
      ++it1;
    } else {
      ++it2;
    }
  }
  return result;
}


double ComputeDistance::calcSummation(const Document& doc1) {

  double result = 0;
  int it1 = 0;

  while (it1 < doc1.iv.size() ) {
    	result += doc1.iv[it1].value;
      ++it1;
    }
  return result;
}

}  // namespace learning_psc

//using namespace learning_psc;
//using namespace std;
//int _main_(int argc, char** argv) {
//  MPI_Init(&argc, &argv);
//  int FLAGS_t_nearest_neighbor = 10;
//  string FLAGS_input = "";
//  string FLAGS_output = "";
//  string distance = "";
//
//
//  for (int i = 1; i < argc; ++i) {
//    if (0 == strcmp(argv[i], "--t_nearest_neighbor")) {
//      std::istringstream(argv[i+1]) >> FLAGS_t_nearest_neighbor;
//      ++i;
//    } else if (0 == strcmp(argv[i], "--input")) {
//      FLAGS_input = argv[i+1];
//      ++i;
//    } else if (0 == strcmp(argv[i], "--output")) {
//      FLAGS_output = argv[i+1];
//      ++i;
//    }
//    else if (0 == strcmp(argv[i], "--dist")) {
//    	distance = argv[i+1];
//          ++i;
//       }
//  }
//
//  if (FLAGS_input == "") {
//    cerr << "--input must not be empty" << endl;
//    MPI_Finalize();
//    return 1;
//  }
//  if (FLAGS_output == "") {
//    cerr << "--output must not be empty" << endl;
//    MPI_Finalize();
//    return 1;
//  }
//
//  ComputeDistance compute(FLAGS_t_nearest_neighbor);
//  compute.Read(FLAGS_input);
//  if (0 == strcmp(distance.c_str(), "e")) {
//		// do euclidean distance
//
//	  compute.ParallelComputeEuclidean();
//	} else if (0 == strcmp(distance.c_str(), "b")) {
//		// do city block distance
//		compute.ParallelComputeManhattan();
//
//	} else if (0 == strcmp(distance.c_str(), "c")) {
//		// do correlation distance
//		compute.ParallelComputeCorrelation();
//
//	} else if (0 == strcmp(distance.c_str(), "o")) {
//		// do cosine distance
//		compute.ParallelComputeCosine();
//
//	}else {
//		// do euclidean distance by default
//		compute.ParallelComputeEuclidean();
//	}
//
//
//  compute.Write(FLAGS_output);
//  MPI_Finalize();
//  return 0;
//}
