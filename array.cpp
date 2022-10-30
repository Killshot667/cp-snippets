//  checki f an array is a permutation

bool permutation(int arr[], int n) 
{ 
    // Set to check the count 
    // of non-repeating elements 
    set<int> hash; 
  
    int maxEle = 0; 
  
    for (int i = 0; i < n; i++) { 
  
        // Insert all elements in the set 
        hash.insert(arr[i]); 
  
        // Calculating the max element 
        maxEle = max(maxEle, arr[i]); 
    } 
  
    if (maxEle != n) 
        return false; 
  
    // Check if set size is equal to n 
    if (hash.size() == n) 
        return true; 
  
    return false; 
} 

// inversion count using pbds (adjacent swaps)
int getInvCount(int arr[], int n)
{
    int key;
    // Initialise the ordered_set
    ordered_set set1;
 
    // Insert the first
    // element in set
    set1.insert(arr[0]);
 
    // Initialise inversion
    // count to zero
    int invcount = 0;
 
    // Finding the inversion
    // count for current element
    for (int i = 1; i < n; i++) {
        set1.insert(arr[i]);
 
        // Number of elements strictly
        // less than arr[i]+1
        key = set1.order_of_key(arr[i] + 1);
 
        // Difference between set size
        // and key will give the
        // inversion count
        invcount += set1.size() - key;
    }
    return invcount;
}