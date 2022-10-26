#include <vector>
#include <iostream>

#include <chrono>

int main()
{
	using std::cout, std::endl;

	std::vector<int> vec;

	cout << "push back s/t ... " << endl;
	vec.push_back(10);
	vec.push_back(10);
	vec.push_back(10);
	vec.push_back(10);
	vec.push_back(10);
	cout << "size check : ";
	cout << vec.size() << endl;

	for( auto x : vec )
	{	cout << x << endl;
	}

	cout << "clear" << endl;
	vec.clear();
	cout << "size check : ";
	cout << vec.size() << endl;

	for( auto x : vec )
	{	cout << x << endl;
	}

	cout << "push back s/t ... " << endl;
	vec.push_back(11);
	vec.push_back(11);
	vec.push_back(11);
	vec.push_back(11);
	cout << "size check : ";
	cout << vec.size() << endl;
	
	for( auto x : vec )
	{	 cout << x << endl;
	}

	cout << "Time Measure" << endl;

	auto start = std::chrono::system_clock::now();
	for(int i=0;i<999999999;i++) {}
	auto end   = std::chrono::system_clock::now();

	std::chrono::duration<double> wtime = end - start;
	cout << wtime.count() << " s\n";


	return 0;
}
