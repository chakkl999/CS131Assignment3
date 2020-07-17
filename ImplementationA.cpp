#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include "mpi.h"

using namespace std;

int* getAdjacencyMatrix(char *fileName, int numtasks);
int* getImage(char *fileName, int &image_size, int rank);
int* computeHistogram(int *partialImage, int size);
void addHistogram(int *first, int *second);
void outputHistogram(const char *fileName, int *histogram);

int main(int argc, char *argv[])
{
    int rank, numtasks;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int image_size;

    int *fullMatrix;
    int *image;
    int *scatterImageSize;

    int *partialImage;
    int *matrix;

    if(rank == 0)
    {
        fullMatrix = getAdjacencyMatrix(argv[2], numtasks);
        image = getImage(argv[1], image_size, rank);
        if(image == nullptr)
            MPI_Finalize();
        scatterImageSize = new int[numtasks];
        for(int i = 0; i < numtasks; ++i)
            scatterImageSize[i] = image_size;
    }
    //if(rank == 0)
    //  {
    //	cout << "Full matrix: " << endl;
    //	for(int i = 0; i < numtasks*numtasks; i+=numtasks)
    //	{
    //	    cout << fullMatrix[i] << " " << fullMatrix[i+1] << " " << fullMatrix[i+2] << " " << fullMatrix[i+3] << endl;
    //	}
    //  }

    MPI_Scatter(scatterImageSize, 1, MPI_INT, &image_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

    partialImage = new int[image_size/numtasks];
    matrix = new int[numtasks];

    MPI_Scatter(fullMatrix, numtasks, MPI_INT, matrix, numtasks, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(image, image_size/numtasks, MPI_INT, partialImage, image_size/numtasks, MPI_INT, 0, MPI_COMM_WORLD);

    //if(rank == 0)
    //    cout << "Finished scattering." << endl;
    int *partialHistogram = computeHistogram(partialImage, image_size/numtasks);
    //string outf = to_string(rank) + ".txt";
    //outputHistogram(outf.c_str(), partialHistogram);
    int message[256];
    for(int i = 0; i < 256; ++i)
        message[i] = -1;

    //    cout << "Node " << rank << "'s adjacency matrix: " << endl;
    //    for(int i = 0; i < numtasks; ++i)
    //	cout << matrix[i] << " ";
    //    cout << endl;

    int parentRank = -1, index = 0, lastVisited = -1;
    if(rank != 2)
    {
        //cout << "Node " << rank << " is waiting for a parent." << endl;
        MPI_Recv(&message, 256, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
        parentRank = status.MPI_SOURCE;
        //cout << "Node " << rank << "'s parent is " << parentRank << endl;
    }
    while(index < numtasks)
    {
        if(matrix[index] == 1 && index > lastVisited && index != parentRank)
        {
            //cout << "Node " << rank << " sending message to " << index << endl;
            MPI_Send(message, 256, MPI_INT, index, 0, MPI_COMM_WORLD);
            lastVisited = index;
            MPI_Recv(message, 256, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
            //cout << "Node " << rank << " recieved message from " << status.MPI_SOURCE << endl;
            //printHistogram(message);
            if(message[0] != -1)
            {
                //cout << "Node " << rank << " is adding histogram from " << status.MPI_SOURCE << endl;
                addHistogram(partialHistogram, message);
            }
            index = 0;
        }
        else
            index++;
    }
    if(rank != 2)
    {
        //cout << "Node " << rank << " send histogram back to parent " << parentRank << endl;
        MPI_Send(partialHistogram, 256, MPI_INT, parentRank, 0, MPI_COMM_WORLD);
    }
    if(rank == 2)
        outputHistogram(argv[3], partialHistogram);
    MPI_Finalize();
    return 0;
}

int* getAdjacencyMatrix(char *fileName, int numtasks)
{
    ifstream in(fileName);
    string line;
    vector<int> matrix;
    int temp;
    while(getline(in, line))
    {
        //	cout << "Line: " << line << endl;
        stringstream ss;
        ss << line;
        //	cout << "Stringstream: " << ss.str() << endl;
        for(int i = 0; i < numtasks; ++i)
        {
            ss >> temp;
            // cout << temp << " ";
            matrix.emplace_back(temp);
        }
        //cout << endl;
    }
    //cout << "Matrix size: " << matrix.size() << endl;
    //for(int i = 0; i < matrix.size(); i+=4)
    //    cout << matrix[i] << " " << matrix[i+1] << " " << matrix[i+2] << " " << matrix[i+3] << endl;
    int* arr = new int[matrix.size()];
    std::copy(matrix.begin(), matrix.end(), arr);
    in.close();
    return arr;
}

int* getImage(char *fileName, int &image_size, int rank)
{
    std::ifstream file(fileName);
    std::string workString;
    int image_width, image_height, image_maxShades;
    /* Remove comments '#' and check image format */
    while(std::getline(file,workString))
    {
        if( workString.at(0) != '#' ){
            if( workString.at(1) != '2' ){
                std::cout << "Input image is not a valid PGM image" << std::endl;
                return nullptr;
            } else {
                break;
            }
        } else {
            continue;
        }
    }
    /* Check image size */
    while(std::getline(file,workString))
    {
        if( workString.at(0) != '#' ){
            std::stringstream stream(workString);
            int n;
            stream >> n;
            image_width = n;
            stream >> n;
            image_height = n;
            break;
        } else {
            continue;
        }
    }

    /* Check image max shades */
    while(std::getline(file,workString))
    {
        if( workString.at(0) != '#' ){
            std::stringstream stream(workString);
            stream >> image_maxShades;
            break;
        } else {
            continue;
        }
    }
    /* Fill input image matrix */
    int pixel_val, index = 0;
    image_size = image_height*image_width;
    int *inputImage = new int[image_size];
    for( int i = 0; i < image_height; i++ )
    {
        if( std::getline(file,workString) && workString.at(0) != '#' ){
            std::stringstream stream(workString);
            for( int j = 0; j < image_width; j++ ){
                if( !stream )
                    break;
                stream >> pixel_val;
                inputImage[index++] = pixel_val;
            }
        } else {
            continue;
        }
    }
    file.close();
    return inputImage;
}

int* computeHistogram(int *partialImage, int size)
{
    int *histogram;
    histogram = new int[256];
    for(int i = 0; i < 256; ++i)
        histogram[i] = 0;
    for(int i = 0 ; i < size; ++i)
    {
        if(partialImage[i] < 0 || partialImage[i] > 255)
            continue;
        histogram[partialImage[i]]++;
    }
    return histogram;
}

void addHistogram(int *first, int *second)
{
    for(int i = 0; i < 256; ++i)
    {
        first[i] += second[i];
        second[i] = -1;
    }
}

void outputHistogram(const char *fileName, int *histogram)
{
    std::ofstream outfile;
    outfile.open(fileName);
    for(int i = 0; i < 255; ++i)
        outfile << histogram[i] << std::endl;
    outfile << histogram[255];
    outfile.close();
}
