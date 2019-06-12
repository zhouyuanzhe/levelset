#include <iostream>
#include <string>
#include <fstream>


using namespace std;

bool deb=1;
void readHBF( std::string const& s )
{
 
    ifstream in( s, ios::binary );
      
    if ( !in )
    {
        cout << "Error opening file " << s.c_str() << endl;
        exit(0);
    }

    string header;
    char ch;
    size_t count = 0;
    while ((ch = in.get()) != '\0')
    {
        header += ch;
        ++count;
    }

    int32_t version = 0;
    in.read( (char*)&version, sizeof( int32_t ) );
    int32_t rows=0, cols=0;
    in.read( (char*)&rows, sizeof( int32_t ) );
    in.read( (char*)&cols, sizeof( int32_t ) );
    size_t size=rows*cols;
    if ( deb)
        cout << "rows: " << rows << " , cols: " << cols << " , size: " << size << endl;
	float *x=new float[size];
    in.read( (char*)x, size*sizeof(float) );
    string nom;
    int i=0;
    while (s[i]!='.') 
	{
		nom+=s[i];
		++i;
	}
	nom+=".txt";
    ofstream out(nom,ios::out | ios::trunc);
    out<<version<<' ';
    out<<rows<<' ';
    out<<cols<<' ';
    for (int i=0;i<size;++i)
	{
		out<<x[i]<<' ';
	}
   	delete x;
   	in.close();
   	out.close();
    return;
}
   

int main()
{
	bool ok;
	string s;
	do
	{
	cout<<"Give the name of the file to convert (without '.hbf')"<<endl;
	cin>>s;
	cout<<endl;
	s+=".hbf";
	readHBF( s );
	cout<<"Do you want to convert an other one (0=no, 1=yes) ?"<<endl;
	cin>>ok;
	cout<<endl;
    }while(ok);
    cout<<"Thank you for using this application."<<endl;
	return 0;
}




