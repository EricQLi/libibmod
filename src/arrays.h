/*
 * arrays.h
 *
 *  Created on: 1 maj 2014
 *      Author: s12kb2
 */

#ifndef ARRAYS_H_
#define ARRAYS_H_

#include <cstdio>
#include <cstdlib>


/**
 * Resizes a C array
 * @param arr a pointer
 * @param new_size
 * @note  'arr' must be either uninitialized or initialized with 'malloc'
 * @return TRUE if resizing succeeded
 */
template <class T>
bool arr_resize(T* & arr, const size_t new_size) {
	T* tmp = static_cast<T*> (realloc(arr, (new_size) * sizeof(T)));
	if (tmp != NULL) {
		arr = tmp;
		return true;
	} else {
		//free(arr);
		fprintf(stderr, "Error: in 'arr_resize', realloc failed (new_size = %lu) \n", (unsigned long) new_size);
		return false;
	}
}

template <class T>
bool arr_resize2(T* & arr, size_t n, size_t & size) {
	if(n >= size) {
		size_t NEW_size = size * 2;
		if(arr_resize(arr, NEW_size)) {
			size = NEW_size;
		} else {
			fprintf(stderr, "Error: in 'arr_resize2', cannot resize array to %lu \n", (unsigned long) NEW_size);
			return false;
		}
	}
	return true;
}

template <class T>
bool arr_push(T* & arr, const T & item, size_t & pos, size_t & size) {
	if(pos >= size) {
		size_t NEW_size = (size * 2) + 1;
		if(arr_resize(arr, NEW_size)) {
			size = NEW_size;
		} else {
			fprintf(stderr, "Error: in 'arr_push', cannot resize array to %lu \n", (unsigned long) NEW_size);
			return false;
		}
	}
	arr[pos++] = item;
	return true;
}


template<class T>
class ExtArray {
private:
	size_t _size;
	size_t _pos;
	T* _iterator;
public:
	T *values;

	ExtArray(size_t sz) :  _size(sz),  _pos(0), _iterator(NULL), values(NULL) {
		if(_size != 0) arr_resize(values, _size);
	}

	ExtArray(void) :  _size(0), _pos(0), _iterator(NULL), values(NULL) {
	}

	void push_back(const T & element) {
		arr_push(values, element, _pos, _size);
	}

	size_t length(void) const {
		return _pos;// XXX: new items are inserted at @_pos. Same as array length.
	}

	T & operator[](const size_t idx) {
//		printf("operator[] #1 \n");
		while(idx >= _pos) {
//			printf("pushing one element to %d \n", _pos);
			arr_push(values, static_cast<T>(0), _pos, _size);
		}
		return values[idx];
	}


    const T & operator[](const size_t idx) const {
    	printf("using const operator[] #2 \n");
		while(idx >= _pos) {
//			printf("pushing one element to %lu \n", _pos);
			arr_push(values, static_cast<T>(0), _pos, _size);
		}
        return const_cast<T&>(values[idx]);
    };

	T & at(const size_t pos)  {
	   return values[pos];
	}

	void clear(void) {
		for (size_t i = 0; i < _pos; ++i) {
			values[i] = static_cast<T>(0);
		}
		_pos = 0;
	}

	size_t pop_back(void){
		values[--_pos] = static_cast<T>(0);
		return _pos;
	}

	T & remove_at(const size_t pos) {
		if (pos >= _pos) {
			errno = 10001;
			perror("pos > array_size");
			exit(10001);
		}
		T & el = values[pos];
		values[pos] = values[--_pos];
		return el;
	}

	bool remove(const T & element) {
		for (size_t i = 0; i < _pos; ++i) {
			if (values[i] == element) {
				remove_at(i);
				return true;
			}
		}
		return false;
	}

	~ExtArray() {
		free(values);
	}


	T* reset(void) {
		_iterator = values;
		return _iterator;
	}

	T* next(void)  {
		if(_iterator >= values + _pos - 1) return NULL;
		return ++_iterator;
	}


};


#endif /* ARRAYS_H_ */
