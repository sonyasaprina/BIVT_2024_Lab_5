using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Numerics;
using System.Reflection;
using System.Reflection.Metadata;
using System.Text;
using static System.Runtime.InteropServices.JavaScript.JSType;

public class Program
{
    public static void Main()
    {
        Program program = new Program();
    }
    #region Level 1
    public long Task_1_1(int n, int k)
    {
        long answer = 0;

        // code here

        if (k == 0 || k > 0 && k == n) answer = 1;
        else if (k > 0 && k < n) answer = Combinations(n, k);
        else answer = 0;
        // end
        return answer;

    }
    static int Combinations(int n, int k)
    {
        return Factorial(n) / (Factorial(k) * (Factorial(n - k)));
    }

    static int Factorial(int n)
    {
        int f = 1;
        for (int i = 2; i <= n; i++)
        {
            f *= i;
        }
        return f;
    }
    public int Task_1_2(double[] first, double[] second)
    {
        int answer = 0;

        // code here

        double a = first[0], a1 = first[1], a2 = first[2], b = second[0], b1 = second[1], b2 = second[2];
        if (a + a1 <= a2 || a + a2 <= a1 || a2 + a1 <= a || b + b1 <= b2 || b1 + b2 <= b || b2 + b <= b1) return -1;
        if (GeronArea(a, a1, a2) > GeronArea(b, b1, b2)) return 1;
        else if (GeronArea(a, a1, a2) < GeronArea(b, b1, b2)) return 2;
        else if (GeronArea(a, a1, a2) == GeronArea(b, b1, b2)) return 0;
        // end

        // first = 1, second = 2, equal = 0, error = -1
        return answer;
    }
    static double GeronArea(double a, double b, double c)
    {
        double p = (a + b + c) / 2;
        return Math.Sqrt(p * (p - a) * (p - b) * (p - c));
    }

    public int Task_1_3a(double v1, double a1, double v2, double a2, int time)
    {
        int answer = 0;

        // code here

        if (GetDistance(v1, a1, time) > GetDistance(v2, a2, time)) return 1;
        if (GetDistance(v1, a1, time) < GetDistance(v2, a2, time)) return 2;
        else if (GetDistance(v1, a1, time) >= GetDistance(v2, a2, time)) return 0;

        // end

        // first = 1, second = 2, equal = 0
        return answer;
    }

    public int Task_1_3b(double v1, double a1, double v2, double a2)
    {
        int answer = 0;

        // code here
        for (int t = 1; ; t++)
        {
            if (GetDistance(v1, a1, t) <= GetDistance(v2, a2, t)) return t;
        }
        // end

        return answer;
    }
    static double GetDistance(double v, double a, int t)
    {
        return v * t + a * t * t / 2;
    }

    #endregion

    #region Level 2
    public void Task_2_1(int[,] A, int[,] B)
    {
        // code here

        // create and use FindMaxIndex(matrix, out row, out column);

        // end
    }

    public void Task_2_2(double[] A, double[] B)
    {
        // code here
        int ai = 0, bi = 0;
        FindMaxIndex(A, out ai);
        FindMaxIndex(B, out bi);
        if ((A.Length - ai) > (B.Length - bi))
        {
            int k = 0;
            double s = 0, sr;
            for (int i = ai + 1; i < A.Length; i++)
            {
                k++;
                s += A[i];
            }
            sr = s / k;
            A[ai] = sr;
        }
        else
        {
            int k = 0;
            double s = 0, sr;
            for (int i = bi + 1; i < B.Length; i++)
            {
                k++;
                s += B[i];
            }
            sr = s / k;
            B[bi] = sr;


            // end
        }
    }
    static void FindMaxIndex(double[] array, out int imax)
    {
        imax = 0;
        double max = array[0];
        for (int i = 0; i < array.Length; i++)
        {
            if (array[i] > max)
            {
                max = array[i];
                imax = i;
            }
        }
    }

    public void Task_2_3(ref int[,] B, ref int[,] C)
    {
        // code here

        //  create and use method FindDiagonalMaxIndex(matrix);

        // end
    }

    public void Task_2_4(int[,] A, int[,] B)
    {
        // code here
        int n = 0, k = 0, ai = 0, bj = 0;
        FindDiagonalMaxIndex(A, out ai);
        FindDiagonalMaxIndex(B, out bj);
        while (true)
        {
            if (k == 5) break;
            n = A[ai, k];
            A[ai, k] = B[k, bj];
            B[k, bj] = n;
            k++;
        }

        // end
    }
    static void FindDiagonalMaxIndex(int[,] matrix, out int imax)
    {
        int max = matrix[0, 0];
        imax = 0;
        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            if (matrix[i, i] > max)
            {
                max = matrix[i, i];
                imax = i;
            }
        }
    }

    public void Task_2_5(int[,] A, int[,] B)
    {
        // code here

        // create and use FindMaxInColumn(matrix, columnIndex, out rowIndex);

        // end
    }

    public void Task_2_6(ref int[] A, int[] B)
    {
        // code here
        int[] a = new int[A.Length + B.Length - 2];
        int j = 0;
        for (int i = 0; i < a.Length; i++)
        {
            if (i < A.Length - 1)
            {
                a[i] = (DeleteElement(A, FindMax(A)))[i];
            }
            else
            {
                a[i] = (DeleteElement(B, FindMax(B)))[j];
                j++;
            }
        }
        A = a;
        // end
    }
    static int FindMax(int[] array)
    {
        int imax = 0;
        int max = array[0];
        for (int i = 0; i < array.Length; i++)
        {
            if (array[i] > max)
            {
                max = array[i];
                imax = i;
            }
        }
        return imax;
    }
    static int[] DeleteElement(int[] array, int index)
    {
        int[] x = new int[array.Length - 1];
        for (int i = 0; i < x.Length; i++)
        {
            if (i < index)
            {
                x[i] = array[i];
            }
            else
            {
                x[i] = array[i + 1];
            }
        }
        return x;
    }

    public void Task_2_7(ref int[,] B, int[,] C)
    {
        // code here

        // create and use CountRowPositive(matrix, rowIndex);
        // create and use CountColumnPositive(matrix, colIndex);

        // end
    }

    public void Task_2_8(int[] A, int[] B)
    {
        // code here
        A = SortArrayPart(A, FindMax(A));
        B = SortArrayPart(B, FindMax(B));
        // end
    }
    static int[] SortArrayPart(int[] array, int startIndex)
    {
        for (int i = startIndex; i < array.Length; i++)
        {
            int key = array[i], j = i - 1;
            while (j >= startIndex + 1 && array[j] > key)
            {
                array[j + 1] = array[j];
                j--;
            }
            array[j + 1] = key;
        }
        return array;
    }

    public int[] Task_2_9(int[,] A, int[,] C)
    {
        int[] answer = default(int[]);

        // code here

        // create and use SumPositiveElementsInColumns(matrix);

        // end

        return answer;
    }

    public void Task_2_10(ref int[,] matrix)
    {
        // code here

        int min = 100000, max = -100000, jmax = 0, jmin = 0;
        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            for (int j = i + 1; j < matrix.GetLength(1); j++)
            {
                if (min > matrix[i, j])
                {
                    min = matrix[i, j];
                    jmin = j;
                }
            }
            for (int j = 0; j <= i; j++)
            {
                if (matrix[i, j] > max)
                {
                    max = matrix[i, j];
                    jmax = j;
                }
            }
        }
        if (jmin > jmax)
        {
            matrix = RemoveColumn(matrix, jmin);
            matrix = RemoveColumn(matrix, jmax);
        }
        else if (jmax > jmin)
        {
            matrix = RemoveColumn(matrix, jmax);
            matrix = RemoveColumn(matrix, jmin);
        }
        else matrix = RemoveColumn(matrix, jmax);
        // end
    }
    static int[,] RemoveColumn(int[,] matrix, int columnIndex)
    {
        int[,] a = new int[matrix.GetLength(0), matrix.GetLength(1) - 1];
        for (int j = 0; j < matrix.GetLength(1) - 1; j++)
        {
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                if (j < columnIndex) a[i, j] = matrix[i, j];
                else a[i, j] = matrix[i, j + 1];
            }
        }
        return a;
    }

    public void Task_2_11(int[,] A, int[,] B)
    {
        // code here

        // use FindMaxIndex(matrix, out row, out column); from Task_2_1

        // end
    }
    public void Task_2_12(int[,] A, int[,] B)
    {
        // code here

        int m = FindMaxColumnIndex(A);
        int n = FindMaxColumnIndex(B);
        for (int i = 0; i < A.GetLength(0); i++)
        {
            int t = A[i, m];
            A[i, m] = B[i, n];
            B[i, n] = t;
        }
        // end
    }
    static int FindMaxColumnIndex(int[,] matrix)
    {
        int jmax = 0, max = -1000000; ;
        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            for (int j = 0; j < matrix.GetLength(1); j++)
            {
                if (matrix[i, j] > max)
                {
                    max = matrix[i, j];
                    jmax = j;
                }
            }
        }
        return jmax;
    }
        public void Task_2_13(ref int[,] matrix)
    {
        // code here

        // create and use RemoveRow(matrix, rowIndex);

        // end
    }

    public void Task_2_14(int[,] matrix)
    {
        // code here
        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            SortRow(matrix, i);
        }
        // end
    }
    static int[,] SortRow(int[,] matrix, int rowIndex)
    {
        for (int i = 0; i < matrix.GetLength(1); i++)
        {
            for (int j = 0; j < matrix.GetLength(1) - i - 1; j++)
            {
                if (matrix[rowIndex, j] > matrix[rowIndex, j + 1])
                {
                    int t = matrix[rowIndex, j];
                    matrix[rowIndex, j] = matrix[rowIndex, j + 1];
                    matrix[rowIndex, j + 1] = t;
                }
            }
        }
        return matrix;
    }

    public int Task_2_15(int[,] A, int[,] B, int[,] C)
    {
        int answer = 0;

        // code here

        // create and use GetAverageWithoutMinMax(matrix);

        // end

        // 1 - increasing   0 - no sequence   -1 - decreasing
        return answer;
    }

    public void Task_2_16(int[] A, int[] B)
    {
        // code here

        SortNegative(A);
        SortNegative(B);

        // end
    }
    static int[] SortNegative(int[] array)
    {
        int k = 0;
        int[] a = new int[array.Length];
        for (int i = 0; i < array.Length; i++)
        {
            if (array[i] < 0)
            {
                a[k] = array[i];
                k++;
            }
        }
        for (int i = 0; i < k; i++)
        {
            for (int j = 0; j < k - i - 1; j++)
            {
                if (a[j] < a[j + 1])
                {
                    int t = a[j];
                    a[j] = a[j + 1];
                    a[j + 1] = t;
                }
            }
        }
        int n = 0;
        for (int i = 0; i < array.Length; i++)
        {
            if (array[i] < 0)
            {
                array[i] = a[n];
                n++;
            }
        }
        return array;
    }

    public void Task_2_17(int[,] A, int[,] B)
    {
        // code here

        // create and use SortRowsByMaxElement(matrix);

        // end
    }

    public void Task_2_18(int[,] A, int[,] B)
    {
        // code here

        SortDiagonal(A);
        SortDiagonal(B);

        // end
    }
    static int[,] SortDiagonal(int[,] matrix)
    {
        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            int key = matrix[i, i], j = i - 1;
            while (j >= 0 && matrix[j, j] > key)
            {
                matrix[j + 1, j + 1] = matrix[j, j];
                j--;
            }
            matrix[j + 1, j + 1] = key;
        }
        return matrix;
    }

    public void Task_2_19(ref int[,] matrix)
    {
        // code here

        // use RemoveRow(matrix, rowIndex); from 2_13

        // end
    }
    public void Task_2_20(ref int[,] A, ref int[,] B)
    {
        // code here

        for (int j = A.GetLength(1) - 1; j >= 0; j--)
        {
            int k = 0;
            for (int i = 0; i < A.GetLength(0); i++)
            {
                if (A[i, j] == 0) k++;
            }
            if (k == 0) A = RemoveColumn(A, j);
        }
        for (int j = B.GetLength(1) - 1; j >= 0; j--)
        {
            int k = 0;
            for (int i = 0; i < B.GetLength(0); i++)
            {
                if (B[i, j] == 0) k++;
            }
            if (k == 0) B = RemoveColumn(B, j);
        }


        // end
    }

    public void Task_2_21(int[,] A, int[,] B, out int[] answerA, out int[] answerB)
    {
        answerA = null;
        answerB = null;

        // code here

        // create and use CreateArrayFromMins(matrix);

        // end
    }

    public void Task_2_22(int[,] matrix, out int[] rows, out int[] cols)
    {
        rows = null;
        cols = null;

        // code here
        rows = new int[matrix.GetLength(0)];
        cols = new int[matrix.GetLength(1)];


        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            rows[i] = CountNegativeInRow(matrix, i);
        }

        for (int j = 0; j < matrix.GetLength(1); j++)
        {
            cols[j] = FindMaxNegativePerColumn(matrix, j);
        }
        // end
    }
    static int CountNegativeInRow(int[,] matrix, int rowIndex)
    {
        int k = 0;
        for (int j = 0; j < matrix.GetLength(1); j++)
        {
            if (matrix[rowIndex, j] < 0) k++;
        }
        return k;
    }
    static int FindMaxNegativePerColumn(int[,] matrix, int colIndex)
    {
        int max = -100000;
        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            if (matrix[i, colIndex] > max && matrix[i, colIndex] < 0) max = matrix[i, colIndex];
        }
        return max;
    }

    public void Task_2_23(double[,] A, double[,] B)
    {
        // code here

        // create and use MatrixValuesChange(matrix);

        // end
    }

    public void Task_2_24(int[,] A, int[,] B)
    {
        // code here
        int i = 0, j = 0;
        FindMaxIndex(A, out i, out j);
        A = SwapColumnDiagonal(A, j);
        FindMaxIndex(B, out i, out j);
        B = SwapColumnDiagonal(B, j);
        // end
    }
    static int FindMaxIndex(int[,] matrix, out int row, out int column)
    {
        row = 0; column = 0;
        int max = -100000;
        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            for (int j = 0; j < matrix.GetLength(1); j++)
            {
                if (matrix[i, j] > max)
                {
                    max = matrix[i, j];
                    row = i;
                    column = j;
                }
            }
        }
        return column;
    }
    static int[,] SwapColumnDiagonal(int[,] matrix, int columnIndex)
    {

        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            int t = matrix[i, i];
            matrix[i, i] = matrix[i, columnIndex];
            matrix[i, columnIndex] = t;
        }
        return matrix;

    }

    public void Task_2_25(int[,] A, int[,] B, out int maxA, out int maxB)
    {
        maxA = 0;
        maxB = 0;

        // code here

        // create and use FindRowWithMaxNegativeCount(matrix);
        // in FindRowWithMaxNegativeCount create and use CountNegativeInRow(matrix, rowIndex); like in 2_22

        // end
    }

    public void Task_2_26(int[,] A, int[,] B)
    {
        // code here

        int imaxA = FindRowWithMaxNegativeCount(A);
        int imaxB = FindRowWithMaxNegativeCount(B);
        for (int i = 0; i < A.GetLength(0); i++)
        {
            for (int j = 0; j < A.GetLength(1); j++)
            {
                int t = A[imaxA, j];
                A[imaxA, j] = B[imaxB, j];
                B[imaxB, j] = t;
            }
        }
        // end
    }
    static int FindRowWithMaxNegativeCount(int[,] matrix)
    {
        int max = 0, imax = 0;
        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            if (max < CountNegativeInRow(matrix, i))
            {
                max = CountNegativeInRow(matrix, i);
                imax = i;
            }
        }
        return imax;
    }

    public void Task_2_27(int[,] A, int[,] B)
    {
        // code here

        // create and use FindRowMaxIndex(matrix, rowIndex, out columnIndex);
        // create and use ReplaceMaxElementOdd(matrix, row, column);
        // create and use ReplaceMaxElementEven(matrix, row, column);

        // end
    }
    static int FindSequence(int[] array, int A, int B)
    {
        bool increasing = true, decreasing = true;

        for (int i = A; i < B; i++)
        {
            if (array[i] < array[i + 1]) decreasing = false;

            else if (array[i] > array[i + 1]) increasing = false;
        }
        if (increasing) return 1;
        if (decreasing) return -1;
        else return 0;
    }

    static int[,] FindSequences_(int[] array, int A, int B)
    {
        int k = 0, ind = 0;
        for (int i = A; i < B; i++)
        {
            for (int j = i + 1; j <= B; j++)
            {
                if (FindSequence(array, i, j) != 0) k++;
            }
        }
        int[,] answer = new int[k, 2];
        for (int i = A; i < B; i++)
        {
            for (int j = i + 1; j <= B; j++)
            {
                if (FindSequence(array, i, j) != 0)
                {
                    answer[ind, 0] = i;
                    answer[ind, 1] = j;
                    ind++;
                }
            }
        }
        return answer;
    }
    static int[] FindMaxSequences__(int[,] array)
    {
        int maxi = -10000, start = -1, end = -1;
        for (int i = 0; i < array.GetLength(0); i++)
        {
            int st = array[i, 0], en = array[i, 1], m = en - st;

            if (m > maxi)
            {
                maxi = m;
                start = st;
                end = en;
            }
        }
        int[] max = new int[] { start, end };

        return max;
    }
    public void Task_2_28a(int[] first, int[] second, ref int answerFirst, ref int answerSecond)
    {
        // code here

        answerFirst = FindSequence(first, 0, first.Length - 1);
        answerSecond = FindSequence(second, 0, second.Length - 1);


        // end
    }

    public void Task_2_28b(int[] first, int[] second, ref int[,] answerFirst, ref int[,] answerSecond)
    {
        // code here
        answerFirst = FindSequences_(first, 0, first.Length - 1);
        answerSecond = FindSequences_(second, 0, second.Length - 1);
        // end
    }

    public void Task_2_28c(int[] first, int[] second, ref int[] answerFirst, ref int[] answerSecond)
    {
        // code here

        int[,] a = FindSequences_(first, 0, first.Length - 1);
        answerFirst = FindMaxSequences__(a);
        int[,] b = FindSequences_(second, 0, second.Length - 1);
        answerSecond = FindMaxSequences__(b);

        // end
    }
    #endregion

    #region Level 3
    public void Task_3_1(ref double[,] firstSumAndY, ref double[,] secondSumAndY)
    {
        // code here

        // create and use public delegate SumFunction(x) and public delegate YFunction(x)
        // create and use method GetSumAndY(sFunction, yFunction, a, b, h);
        // create and use 2 methods for both functions calculating at specific x

        // end
    }
    public delegate int[,] SortRowStyle(int[,] matrix, int rowIndex);

    static int[,] SortAscending(int[,] matrix, int rowIndex)
    {
        for (int i = 0; i < matrix.GetLength(1); i++)
        {
            int key = matrix[rowIndex, i], j = i - 1;
            while (j >= 0 && matrix[rowIndex, j] > key)
            {
                matrix[rowIndex, j + 1] = matrix[rowIndex, j];
                j--;
            }
            matrix[rowIndex, j + 1] = key;
        }
        return matrix;
    }

    static int[,] SortDescending(int[,] matrix, int rowIndex)
    {
        for (int i = 0; i < matrix.GetLength(1); i++)
        {
            int key = matrix[rowIndex, i], j = i - 1;
            while (j >= 0 && matrix[rowIndex, j] < key)
            {
                matrix[rowIndex, j + 1] = matrix[rowIndex, j];
                j--;
            }
            matrix[rowIndex, j + 1] = key;
        }
        return matrix;
    }

    public void Task_3_2(int[,] matrix)
    {
        SortRowStyle sortStyle = default(SortRowStyle);

        // code here

        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            if (i % 2 == 0) sortStyle = SortAscending;
            else sortStyle = SortDescending;
            sortStyle(matrix, i);
        }


        // end
    }

    public double Task_3_3(double[] array)
    {
        double answer = 0;
        // SwapDirection swapper = default(SwapDirection); - uncomment me

        // code here

        // create and use public delegate SwapDirection(array);
        // create and use methods SwapRight(array) and SwapLeft(array)
        // create and use method GetSum(array, start, h) that sum elements with a negative indexes
        // change method in variable swapper in the if/else and than use swapper(matrix)

        // end

        return answer;
    }
    public delegate int[] GetTriangle(int[,] matrix);

    static int[] GetUpperTriange(int[,] matrix)
    {
        int k = 0;
        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            for (int j = i; j < matrix.GetLength(1); j++)
            {
                k++;
            }
        }
        int[] array = new int[k];
        int ind = 0;
        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            for (int j = i; j < matrix.GetLength(1); j++)
            {
                array[ind] = matrix[i, j] * matrix[i, j];
                ind++;
            }
        }
        return array;
    }

    static int[] GetLowerTriange(int[,] matrix)
    {
        int k = 0;
        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            for (int j = 0; j <= i; j++)
            {
                k++;
            }
        }
        int[] array = new int[k];
        int ind = 0;
        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            for (int j = 0; j <= i; j++)
            {
                array[ind] = matrix[i, j] * matrix[i, j];
                ind++;
            }
        }
        return array;
    }
    static int GetSum(GetTriangle a, int[,] matrix)
    {
        int sum = 0;
        int[] array = a(matrix);
        for (int i = 0; i < array.Length; i++)
        {
            sum += array[i];
        }
        return sum;
    }

    public int Task_3_4(int[,] matrix, bool isUpperTriangle)
    {
        int answer = 0;

        // code here
        GetTriangle a;
        if (isUpperTriangle) a = GetUpperTriange;
        else a = GetLowerTriange;
        answer = GetSum(a, matrix);

        // end

        return answer;
    }

    public void Task_3_5(out int func1, out int func2)
    {
        func1 = 0;
        func2 = 0;

        // code here

        // use public delegate YFunction(x, a, b, h) from Task_3_1
        // create and use method CountSignFlips(YFunction, a, b, h);
        // create and use 2 methods for both functions

        // end
    }
    public delegate int FindElementDelegate(int[,] matrix);

    static int FindDiagonalMaxIndex(int[,] matrix)
    {
        int max = -100000, d = 0;
        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            if (matrix[i, i] > max)
            {
                max = matrix[i, i];
                d = i;
            }
        }
        return d;
    }
    static int FindFirstRowMaxIndex(int[,] matrix)
    {
        int max = -100000, f = 0;
        for (int i = 0; i < matrix.GetLength(1); i++)
        {
            if (matrix[0, i] > max)
            {
                max = matrix[0, i];
                f = i;
            }
        }
        return f;
    }
    static int[,] SwapColumns(int[,] matrix, FindElementDelegate first, FindElementDelegate second)
    {
        int j = first(matrix), ind = second(matrix);
        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            int t = matrix[i, j];
            matrix[i, j] = matrix[i, ind];
            matrix[i, ind] = t;
        }
        return matrix;
    }

    public void Task_3_6(int[,] matrix)
    {
        // code here

        FindElementDelegate f = FindDiagonalMaxIndex;
        FindElementDelegate s = FindFirstRowMaxIndex;
        matrix = SwapColumns(matrix, f, s);

        // end
    }

    public void Task_3_7(ref int[,] B, int[,] C)
    {
        // code here

        // create and use public delegate CountPositive(matrix, index);
        // use CountRowPositive(matrix, rowIndex) from Task_2_7
        // use CountColumnPositive(matrix, colIndex) from Task_2_7
        // create and use method InsertColumn(matrixB, CountRow, matrixC, CountColumn);

        // end
    }

    public delegate int FindIndex(int[,] matrix);
    static int FindMaxBelowDiagonalIndex(int[,] matrix)
    {
        int max = -10000, ind = 0;
        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            for (int j = 0; j <= i; j++)
            {
                if (matrix[i, j] > max)
                {
                    max = matrix[i, j];
                    ind = j;
                }
            }
        }
        return ind;
    }
    static int FindMinAboveDiagonalIndex(int[,] matrix)
    {
        int min = 10000, ind = 0;
        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            for (int j = i; j < matrix.GetLength(1); j++)
            {
                if (matrix[i, j] < min)
                {
                    min = matrix[i, j];
                    ind = j;
                }
            }
        }
        return ind;
    }
    static int[,] RemoveColumns1(int[,] matrix, FindIndex b, FindIndex a)
    {
        int m = b(matrix), n = a(matrix);
        if (n > m)
        {
            matrix = RemoveColumn(matrix, n);
            matrix = RemoveColumn(matrix, m);
        }
        else if (m > n)
        {
            matrix = RemoveColumn(matrix, m);
            matrix = RemoveColumn(matrix, n);
        }
        else matrix = RemoveColumn(matrix, m);
        return matrix;
    }

    public void Task_3_10(ref int[,] matrix)
    {
        // FindIndex searchArea = default(FindIndex); - uncomment me

        // code here
        FindIndex b = FindMaxBelowDiagonalIndex;
        FindIndex a = FindMinAboveDiagonalIndex;
        matrix = RemoveColumns1(matrix, b, a);

        // end
    }

    public void Task_3_13(ref int[,] matrix)
    {
        // code here

        // use public delegate FindElementDelegate(matrix) from Task_3_6
        // create and use method RemoveRows(matrix, findMax, findMin)

        // end
    }
    public delegate int GetNegativeArray(int[,] matrix, int ind);
    static void FindNegatives(int[,] matrix, GetNegativeArray searcherRows, GetNegativeArray searcherCols, out int[] rows, out int[] cols)
    {
        rows = new int[matrix.GetLength(0)];
        cols = new int[matrix.GetLength(1)];
        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            rows[i] = searcherRows(matrix, i);
        }

        for (int j = 0; j < matrix.GetLength(1); j++)
        {
            cols[j] = searcherCols(matrix, j);
        }
    }
    public void Task_3_22(int[,] matrix, out int[] rows, out int[] cols)
    {

        rows = null;
        cols = null;

        // code here

        GetNegativeArray Rows = CountNegativeInRow;
        GetNegativeArray Cols = FindMaxNegativePerColumn;
        FindNegatives(matrix, Rows, Cols, out rows, out cols);


        // end
    }

    public void Task_3_27(int[,] A, int[,] B)
    {
        // code here

        // create and use public delegate ReplaceMaxElement(matrix, rowIndex, max);
        // use ReplaceMaxElementOdd(matrix) from Task_2_27
        // use ReplaceMaxElementEven(matrix) from Task_2_27
        // create and use method EvenOddRowsTransform(matrix, replaceMaxElementOdd, replaceMaxElementEven);

        // end
    }
    public delegate bool IsSequence(int[] array, int left, int right);
    static bool FindIncreasingSequence(int[] array, int A, int B)
    {
        bool answer = false;
        for (int i = A; i < B; i++)
        {
            if (array[i] < array[i + 1]) answer = true;
            else
            {
                answer = false;
                break;
            }
        }
        return answer;
    }
    static bool FindDecreasingSequence(int[] array, int A, int B)
    {
        bool answer = false;
        for (int i = A; i < B; i++)
        {
            if (array[i] > array[i + 1]) answer = true;
            else
            {
                answer = false;
                break;
            }
        }
        return answer;
    }
    static int DefineSequence(int[] array, IsSequence findIncreasingSequence, IsSequence findDecreasingSequence)
    {
        if (findIncreasingSequence(array, 0, array.Length - 1)) return 1;
        if (findDecreasingSequence(array, 0, array.Length - 1)) return -1;
        else return 0;
    }
    static int[] FindLongestSequence(int[] array, IsSequence sequence)
    {
        int[] a = new int[2];
        int max = -1;
        for (int i = 0; i < array.Length - 1; i++)
        {
            for (int j = i + 1; j < array.Length; j++)
            {
                if (sequence(array, i, j))
                {
                    int b = j - i;
                    if (b > max)
                    {
                        max = b;
                        a[0] = i;
                        a[1] = j;
                    }
                }
            }
        }
        return a;
    }

    public void Task_3_28a(int[] first, int[] second, ref int answerFirst, ref int answerSecond)
    {
        // code here

        IsSequence m = FindIncreasingSequence;
        IsSequence n = FindDecreasingSequence;
        answerFirst = DefineSequence(first, m, n);
        answerSecond = DefineSequence(second, m, n);


        // end
    }

    public void Task_3_28c(int[] first, int[] second, ref int[] answerFirstIncrease, ref int[] answerFirstDecrease, ref int[] answerSecondIncrease, ref int[] answerSecondDecrease)
    {
        // code here

        IsSequence m = FindIncreasingSequence;
        IsSequence n = FindDecreasingSequence;
        answerFirstIncrease = FindLongestSequence(first, m);
        answerFirstDecrease = FindLongestSequence(first, n);
        answerSecondIncrease = FindLongestSequence(second, m);
        answerSecondDecrease = FindLongestSequence(second, n);

        // end
    }
    #endregion
    #region bonus part
    public double[,] Task_4(double[,] matrix, int index)
    {
        // MatrixConverter[] mc = new MatrixConverter[4]; - uncomment me

        // code here

        // create public delegate MatrixConverter(matrix);
        // create and use method ToUpperTriangular(matrix);
        // create and use method ToLowerTriangular(matrix);
        // create and use method ToLeftDiagonal(matrix); - start from the left top angle
        // create and use method ToRightDiagonal(matrix); - start from the right bottom angle

        // end

        return matrix;
    }
    #endregion
}
