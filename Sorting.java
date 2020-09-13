import java.util.Random;

import Plotter.Plotter;


public class Sorting {

	final static int SMALL_DEMO_ARRAY = 10000;
	final static int BIG_DEMO_ARRAY = 100000000;
	final static int BUILDHEAP_VS_QUICK_LENGTH = 15;
	final static int MERGE_VS_QUICK_LENGTH = 16;
	final static int MERGE_VS_QUICK_SORTED_LENGTH = 12;
	final static int HEAP_VS_BUBBLE_LENGTH = 16;
	final static double T = 600.0;
	/**
	 * Sorts a given array using the quick sort algorithm.
	 * At each stage the pivot is chosen to be the righttmost element of the subarray.
	 * 
	 * Should run in average complexity of O(nlog(n)), and worst case complexity of O(n^2)
	 * 
	 * @param arr - the array to be sorted
	 */
	public static void quickSort(double[] arr){
		quickSort(arr, 0, arr.length - 1);
	}

	public static void quickSort(double[] arr, int left, int right){
		if (left < right) {
			int inRightPlace = partition(arr, left, right);
			quickSort(arr, left, inRightPlace - 1);
			quickSort(arr, inRightPlace + 1, right);
		}
	}

	public static int partition(double[] arr, int lead, int right){
		double x = arr[right];
		int tail = lead - 1;
		double temp;
		for (; lead < right; lead++) {
			if (arr[lead] < x) {
				swap (arr, ++tail, lead);
			}
		}
		swap (arr, ++tail, right);
		return tail;
	}

	public static void swap (double[] arr, int a, int b) {
		double temp = arr[b];
		arr[b] = arr[a];
		arr[a] = temp;
	}
	
	/**
	 * Sorts a given array using the merge sort algorithm.
	 * 
	 * Should run in complexity O(nlog(n)) in the worst case.
	 * 
	 * @param arr - the array to be sorted
	 */
	public static void mergeSort(double[] arr){
		double[] temp = new double[arr.length];
		mergeSort(arr, temp, 0, arr.length - 1);

	}

	public static void mergeSort(double[] arr, double[] temp, int left, int right) {
		if (left < right) {
			int mid = (right + left) / 2;
			mergeSort(arr, temp, left, mid);
			mergeSort(arr, temp, mid + 1, right);
			merge(arr, temp, left, mid, right);
		}
	}

	public static void merge (double[] arr, double[] temp, int left, int mid, int right) {
		int rightStart = mid + 1;
		int leftIndex = left;
		int rightIndex = mid + 1;
		int index = left;


		while (left <= mid && rightStart <= right) {
			if (arr[left] < arr[rightStart]) {
				temp[index++] = arr[left++];
			}
			else
				temp[index++] = arr[rightStart++];
		}

		if (left > mid) {
			while (rightStart <= right) {
				temp[index++] = arr[rightStart++];
			}
		}
		else {
			while (left <= mid) {
				temp[index++] = arr[left++];
			}
		}

		for (int i = leftIndex; i < right + 1; i++) {
			arr[i] = temp[i];
		}
	}




	/**
	 * Builds a max-heap from the given array.
	 * That is, when using the resulting array to represent an almost complete binary tree
	 * The value of each node must be larger (or equal) to the value of any of its descendants.
	 * Should run in complexity O(n) in the worst case.
	 * 
	 * @param arr - the given array.
	 */
	public static void buildHeap(double[] arr){
		int n = arr.length;
		for (int i = (n / 2) - 1; i >= 0; i--) {
			percDown(i, arr[i], n, arr);
		}
	}

	public static void percDown(int index, double valOfIndex, int length, double[] arr) {
		if ((2 * index) + 1 >= length) {            //has no children
			return;
		}

		else if ((2 * index) + 2 == length) {        //has one child
			if (valOfIndex < arr[(2 * index) + 1]) {
				arr[index] = arr[(2 * index) + 1];
				arr[(2 * index) + 1] = valOfIndex;
				percDown((2 * index) + 1, valOfIndex, length, arr);
			}
			else {
				return;
			}
		}

		else {                                        //has two children
			if (arr[(2 * index) + 1] < arr[(2 * index) + 2]) {
				if (valOfIndex < arr[(2 * index) + 2]) {
					arr[index] = arr[(2 * index) + 2];
					arr[(2 * index) + 2] = valOfIndex;
					percDown((2 * index) + 2, valOfIndex, length, arr);
				}
			}
			else if (valOfIndex < arr[(2 * index) + 1]) {
				arr[index] = arr[(2 * index) + 1];
				arr[(2 * index) + 1] = valOfIndex;
				percDown((2 * index) + 1, valOfIndex, length, arr);
			}
		}
	}


	/**
	 * Sorts a given array using the heap sort algorithm.
	 * 
	 * Should run in complexity O(nlog(n)) in the worst case.
	 * 
	 * @param arr - the array to be sorted
	 */
	public static void heapSort (double[] arr){
		buildHeap(arr);
		int n = arr.length;
		double x;
		for (int i = n - 1 ; i >= 1; i--) {
			x = arr[0];
			arr[0] = arr[i];
			arr[i] = x;
			percDown(0, arr[0], i, arr);
		}
	}




	/**
	 * Sorts a given array using the bubble sort algorithm.
	 * 
	 * Should run in complexity O(n^2) in the worst case.
	 * 
	 * @param arr - the array to be sorted
	 */
	public static void bubbleSort(double[] arr){
		boolean isSorted = false;
		int length = arr.length - 1;
		while (! isSorted) {
			isSorted = true;
			for (int i = 0; i < length; i++) {
				if (arr[i] > arr[i + 1]) {
					swap(arr, i, i + 1);
					isSorted = false;
				}
			}
			length--;
		}
	}

	public static void main(String[] args) {
		buildHeapVsQuick();
		mergeVsQuick();
		mergeVsQuickOnSortedArray();
		heapSortVsBubble();
	}
	
	/**
	 * Compares the build heap algorithm against quick sort on random arrays
	 */
	public static void buildHeapVsQuick(){
		double[] quickTimes = new double[BUILDHEAP_VS_QUICK_LENGTH];
		double[] heapTimes = new double[BUILDHEAP_VS_QUICK_LENGTH];
		long startTime, endTime;
		Random r = new Random();
		for (int i = 0; i < BUILDHEAP_VS_QUICK_LENGTH; i++) {
			long sumQuick = 0;
			long sumHeap = 0;
			for(int k = 0; k < T; k++){
				int size = (int)Math.pow(2, i);
				double[] a = new double[size];
				double[] b = new double[size];
				for (int j = 0; j < a.length; j++) {
					a[j] = r.nextGaussian() * 5000;
					b[j] = a[j];
				}
				startTime = System.currentTimeMillis();
				quickSort(a);
				endTime = System.currentTimeMillis();
				sumQuick += endTime - startTime;
				startTime = System.currentTimeMillis();
				buildHeap(b);
				endTime = System.currentTimeMillis();
				sumHeap += endTime - startTime;
			}
			quickTimes[i] = sumQuick/T;
			heapTimes[i] = sumHeap/T;
		}
		Plotter.plot("quick sort", quickTimes, "build heap", heapTimes);
	}
	
	/**
	 * Compares the merge sort algorithm against quick sort on random arrays
	 */
	public static void mergeVsQuick(){
		double[] quickTimes = new double[MERGE_VS_QUICK_LENGTH];
		double[] mergeTimes = new double[MERGE_VS_QUICK_LENGTH];
		long startTime, endTime;
		Random r = new Random();
		for (int i = 0; i < MERGE_VS_QUICK_LENGTH; i++) {
			long sumQuick = 0;
			long sumMerge = 0;
			for (int k = 0; k < T; k++) {
				int size = (int)Math.pow(2, i);
				double[] a = new double[size];
				double[] b = new double[size];
				for (int j = 0; j < a.length; j++) {
					a[j] = r.nextGaussian() * 5000;
					b[j] = a[j];
				}
				startTime = System.currentTimeMillis();
				quickSort(a);
				endTime = System.currentTimeMillis();
				sumQuick += endTime - startTime;
				startTime = System.currentTimeMillis();
				mergeSort(b);
				endTime = System.currentTimeMillis();
				sumMerge += endTime - startTime;
			}
			quickTimes[i] = sumQuick/T;
			mergeTimes[i] = sumMerge/T;
		}
		Plotter.plot("quick sort", quickTimes, "merge sort", mergeTimes);
	}

	/**
	 * Compares the merge sort algorithm against quick sort on pre-sorted arrays
	 */
	public static void mergeVsQuickOnSortedArray(){
		double[] quickTimes = new double[MERGE_VS_QUICK_SORTED_LENGTH];
		double[] mergeTimes = new double[MERGE_VS_QUICK_SORTED_LENGTH];
		long startTime, endTime;
		for (int i = 0; i < MERGE_VS_QUICK_SORTED_LENGTH; i++) {
			long sumQuick = 0;
			long sumMerge = 0;
			for (int k = 0; k < T; k++) {
				int size = (int)Math.pow(2, i);
				double[] a = new double[size];
				double[] b = new double[size];
				for (int j = 0; j < a.length; j++) {
					a[j] = j;
					b[j] = j;
				}
				startTime = System.currentTimeMillis();
				quickSort(a);
				endTime = System.currentTimeMillis();
				sumQuick += endTime - startTime;
				startTime = System.currentTimeMillis();
				mergeSort(b);
				endTime = System.currentTimeMillis();
				sumMerge += endTime - startTime;
			}
			quickTimes[i] = sumQuick/T;
			mergeTimes[i] = sumMerge/T;
		}
		Plotter.plot("quick sort on sorted array", quickTimes, "merge sort on sorted array", mergeTimes);
	}

	/**
	 * Compares the bubble sort algorithm against heap sort.
	 */
	public static void heapSortVsBubble(){
		double[] bubbleTimes = new double[HEAP_VS_BUBBLE_LENGTH];
		double[] heapTimes = new double[HEAP_VS_BUBBLE_LENGTH];
		long startTime, endTime;
		Random r = new Random();
		for (int i = 0; i < HEAP_VS_BUBBLE_LENGTH; i++) {
			long sumBubble = 0;
			long sumHeap = 0;
			for (int k = 0; k < T; k++) {
				int size = (int)Math.pow(2, i);
				double[] a = new double[size];
				double[] b = new double[size];
				for (int j = 0; j < a.length; j++) {
					a[j] = r.nextGaussian() * 5000;
					b[j] = a[j];
				}
				startTime = System.currentTimeMillis();
				bubbleSort(a);
				endTime = System.currentTimeMillis();
				sumBubble += endTime - startTime;
				startTime = System.currentTimeMillis();
				heapSort(b);
				endTime = System.currentTimeMillis();
				sumHeap += endTime - startTime;
			}
			bubbleTimes[i] = sumBubble/T;
			heapTimes[i] = sumHeap/T;
		}
		Plotter.plot("bubble sort", bubbleTimes, "heap sort", heapTimes);
	}
}
