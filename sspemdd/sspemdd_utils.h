#ifndef SSPEMDD_UTILS
#define SSPEMDD_UTILS

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>

namespace SSPEMDD_utils{
	// construct all combinations of some parameters
	template< typename T >
	bool next_cartesian(std::vector<T> &vii, std::vector<int> &index_arr, T &cur_vi)
	{
		if (index_arr.size() == 0) { // init
			index_arr.resize(vii.size());
			//for( auto &x : index_arr )
			//	x = 0;
			for (std::vector<int> ::iterator it = index_arr.begin(); it != index_arr.end(); ++it)
				*it = 0;
		}
		if (index_arr[0] == -1)
			return false;
		// get current value
		cur_vi.resize(vii.size());
		for (unsigned i = 0; i < index_arr.size(); ++i)
			cur_vi[i] = vii[i][index_arr[i]];
		// check if last iteration
		bool IsLastValue = true;
		for (unsigned i = 0; i < index_arr.size(); ++i) {
			if (index_arr[i] != vii[i].size() - 1) {
				IsLastValue = false;
				break;
			}
		}
		if (IsLastValue)
			index_arr[0] = -1; // condition of stopping
		else {
			// find last changable row to increase its value
			unsigned last_changable = (unsigned)index_arr.size() - 1;
			while (last_changable != -1){
				if (index_arr[last_changable] < (int)(vii[last_changable].size() - 1))
					break;
				--last_changable;
			}
			index_arr[last_changable]++;
			for (unsigned i = last_changable + 1; i < index_arr.size(); ++i)
				index_arr[i] = 0;
		}

		return true;
	}
};

#endif