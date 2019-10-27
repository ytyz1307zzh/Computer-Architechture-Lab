#include <stdio.h>
int a[6] = {4, 2, 3, 1, 5, 6};
void quickSort(int *num, int left, int right)
{
    if (left >= right)
        return;
    int index = left;
    int pivot = num[left];
    for (int i = left; i <= right; i++) {
        if (num[i] < pivot) {
            int temp = num[i];
            num[i] = num[index];
            num[index] = temp;
            index++;
        }
    }
    num[index] = pivot;
    quickSort(num, left, index - 1);
    quickSort(num, index + 1, right);
}
int main()
{
    quickSort(a, 0, 5);
    return 0;
}
