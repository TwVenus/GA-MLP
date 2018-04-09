import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Random;

/**
 * Created by Venus on 2016/11/19.
 */

public class GADemo {
    public static void main(String args[]) throws IOException {

        //讀取txt資料 ====================================================================================
            int feature = 785; //特增值個數
            int data_num = 60000; //data
            int output_num = 4;
            int ga_num = 10;

            String csvFile = "data.csv";
            String line = "";
            String cvsSplitBy = ",";
            double data_set[][] = new double[data_num][feature];
            int index = 0;
            int clounm = feature;
            try (BufferedReader br = new BufferedReader(new FileReader(csvFile))) {
                while ((line = br.readLine()) != null) {

                    // use comma as separator
                    String[] country = line.split(cvsSplitBy);
                    for (int i = 0; i < clounm; i++) {
                        data_set[index][i] = Double.parseDouble(country[i]);
                    }
                    index++;
                }
            } catch (IOException e) {
                e.printStackTrace();
            }
        //===============================================================================================

        FileWriter fwriter = new FileWriter("output.txt");

        int run_count = 0 ;
        int hidden_num = (int)Math.sqrt(feature*output_num);
        int bytes_num = feature * hidden_num + hidden_num + 1;
        int bits_num = 8 * (feature * hidden_num + hidden_num + 1);

        int [][] binary_weight_new = new int[ga_num][bits_num]; //存新的權重

        while(run_count != 2000) {
            //原始的random基因fitness =========================================================================
            System.out.println(run_count);

            int[][] binary_weight = new int[ga_num][bits_num]; // 10條基因，ga_num個weight，bits_num個bit

            if (run_count == 0) {
                for (int i = 0; i < binary_weight.length; i++) {
                    binary_weight[i] = randomWeight(bits_num);
                }
            } else {
                for (int i = 0; i < binary_weight.length; i++) {
                    for (int j = 0; j < binary_weight[i].length; j++)
                        binary_weight[i][j] = binary_weight_new[i][j];
                }

                //初始化
                for (int i = 0; i < binary_weight_new.length; i++) {
                    for (int j = 0; j < binary_weight_new[i].length; j++)
                        binary_weight_new[i][j] = 0;
                }
            }

            double[][] decimal_weight = new double[ga_num][bytes_num];  //10條基因，bytes_num個十進位的weight
            for (int i = 0; i < decimal_weight.length; i++) {
                decimal_weight[i] = binaryToDecimal(binary_weight[i], bytes_num);
            }

            double[] fitness_original = new double[10];  // 儲存原本10個基因的fitness
            fitness_original = fitnessCompute(data_set, decimal_weight, hidden_num);  //計算fitness

            //看fitness_original
            for (int i = 0; i < fitness_original.length; i++) {
                System.out.print(fitness_original[i] + " ");
            }
            System.out.println();
            //===============================================================================================

            //交配基因fitness =================================================================================
            int[][] binary_weight_mating = new int[ga_num][bits_num]; //儲存交配過後的基因
            binary_weight_mating = mating(binary_weight, bytes_num);

            double[][] decimal_weight_mating = new double[ga_num][bytes_num];  //儲存交配過後的基因, 十進位
            for (int i = 0; i < decimal_weight_mating.length; i++) {
                decimal_weight_mating[i] = binaryToDecimal(binary_weight_mating[i], bytes_num);
            }

            double[] fitness_mating = new double[10];  // 儲存10個交配過基因的fitness
            fitness_mating = fitnessCompute(data_set, decimal_weight_mating, hidden_num);  //計算fitness

            //看fitness_mating
            for (int i = 0; i < fitness_mating.length; i++) {
                System.out.print(fitness_mating[i] + " ");
            }
            System.out.println();
            //===============================================================================================

            //突變基因fitness =================================================================================
            int[][] binary_weight_mutation = new int[10][104]; //儲存突變過後的基因
            binary_weight_mutation = mutation(binary_weight, bytes_num);

            double[][] decimal_weight_mutation = new double[ga_num][bytes_num];  //儲存突變過後的基因, 十進位

            for (int i = 0; i < decimal_weight_mutation.length; i++) {
                decimal_weight_mutation[i] = binaryToDecimal(binary_weight_mutation[i], bytes_num);
            }

            double[] fitness_mutation = new double[10];  // 儲存10個突變過基因的fitness
            fitness_mutation = fitnessCompute(data_set, decimal_weight_mutation, hidden_num);  //計算fitness

            //看fitness_mutation
            for (int i = 0; i < fitness_mutation.length; i++) {
                System.out.print(fitness_mutation[i] + " ");
            }
            System.out.println();
            //=====================================================================================================

            //挑出最好的7條 + 3條隨機 =================================================================================
            double[] fitness_all = new double[30];
            int fitness_all_count = 0;

            for (int i = 0; i < fitness_original.length; i++) {
                fitness_all[fitness_all_count] = fitness_original[i];
                fitness_all_count++;
            }
            for (int i = 0; i < fitness_mating.length; i++) {
                fitness_all[fitness_all_count] = fitness_mating[i];
                fitness_all_count++;
            }
            for (int i = 0; i < fitness_mutation.length; i++) {
                fitness_all[fitness_all_count] = fitness_mutation[i];
                fitness_all_count++;
            }

            int[] max_pointer = new int[10];
            max_pointer = bestOfFitness(fitness_all);  //挑選最好7條，與隨機3條，位置

            for (int i = 0; i < max_pointer.length; i++) {
                System.out.print(max_pointer[i] + " ");
            }
            System.out.println();

            double[] fitness_all_new = new double[ga_num]; //新的基因的fitness
            for (int i = 0; i < max_pointer.length; i++) {
                if (max_pointer[i] >= 0 && max_pointer[i] <= 9) {
                    for (int j = 0; j < bits_num; j++)
                        binary_weight_new[i][j] = binary_weight[max_pointer[i]][j];
                    fitness_all_new[i] = fitness_all[max_pointer[i]];
                }
                if (max_pointer[i] >= 10 && max_pointer[i] <= 19) {
                    for (int j = 0; j < bits_num; j++)
                        binary_weight_new[i][j] = binary_weight_mating[max_pointer[i] - 10][j];
                    fitness_all_new[i] = fitness_all[max_pointer[i]];
                }
                if (max_pointer[i] >= 20 && max_pointer[i] <= 29) {
                    for (int j = 0; j < bits_num; j++)
                        binary_weight_new[i][j] = binary_weight_mutation[max_pointer[i] - 20][j];
                    fitness_all_new[i] = fitness_all[max_pointer[i]];
                }
            }

            double[][] decimal_binary_weight_new = new double[ga_num][bytes_num];  //10條基因，13個十進位的weight
            for (int i = 0; i < decimal_binary_weight_new.length; i++) {
                decimal_binary_weight_new[i] = binaryToDecimal(binary_weight_new[i], bytes_num);
            }

            //====================================================================================================

            //寫檔到txt ========================================================================================
            for (int i = 0; i < decimal_binary_weight_new.length; i++) {
                for(int j = 0; j < decimal_binary_weight_new[i].length; j++){
                    fwriter.write(decimal_binary_weight_new[i][j] + " ");
                }
                fwriter.write(fitness_all_new[i] + "\r\n");
            }

            fwriter.write("\r\n");
            //=====================================================================================================

            run_count++;
        }
        fwriter.close();
    }

    //不可重覆的random function
    public static int [] randomPointer(int amount_of_one,int first,int last){
        Random ran = new Random();
        int pointer [] = new int [amount_of_one];

        // random 1的位置且不能重覆
        int number = 0;
        boolean uni;
        for(int i = 0 ; i < amount_of_one ; i++){
            uni = true;

            while(uni) {
                uni = false;
                number = ran.nextInt(last)+first;
                for(int j = 0; j < pointer.length; j++) {
                    if(pointer[j] == number) {
                        uni = true;
                    }
                }
            }
            pointer[i] = number;
        }
        return pointer;
    }

    //random 基因(binary)
    public static int [] randomWeight(int bits_num){
        Random ran = new Random();
        int random_weight [] = new int [bits_num];
        int amount_of_one = ran.nextInt(bits_num)+1; //random 有幾個1
        int pointer [] = randomPointer(amount_of_one, 1, bits_num);  // 1的位置，random 1~104 不可重覆

        for(int i = 0 ; i < bits_num ; i++){ //初始weight為0
            random_weight[i] = 0;
        }

        //把pointer的部分變成1
        for(int i = 0 ; i < bits_num ; i++){
            for(int j = 0 ; j < pointer.length ; j++){
                if (i == pointer[j]){
                    random_weight[i] = 1;
                }
            }
        }
        return random_weight;
    }

    //把 binary 轉成 decimal
    public static double [] binaryToDecimal(int binaryW[], int bytes_num){
        int [][] weight = new int[bytes_num][8];
        double [] decimal_weight = new double[bytes_num];

        int count = 0 ;
        for(int i = 0 ; i < weight.length ; i++){
            for(int j = 0 ; j < weight[i].length ; j++){
                weight[i][j] = binaryW[count];
                count ++;
            }
        }

        for(int i = 0 ; i < weight.length ; i++){
            double decimal = 0;
            for(int j = 1 ; j < weight[i].length ; j++){
                if(weight[i][j] == 1){
                    decimal += Math.pow(2, -(j-1));
                }
                if(weight[i][0] == 1)
                    decimal *= -1;
            }
            decimal_weight[i] = decimal;
        }

        return decimal_weight;
    }

    //計算fitness
    public static double [] fitnessCompute (double [][] data_set, double [][] decimal_weight, int hidden_num){
        double [] fitness_original = new double[10];
        int bias = -1;
        for(int count_weight = 0 ; count_weight < 10 ; count_weight++) {
            double true_coount = 0;

            for (int i = 0; i < data_set.length; i++) {

                double [] hidden_ans = new double[hidden_num]; //儲存到隱藏層所有結果

                for(int k = 0 ; k < hidden_num ; k++){
                    double num1 = 0;
                    for(int j = 0 ; j < data_set[i].length -1 ; j++){
                        num1 += decimal_weight[count_weight][data_set[i].length*k + j + 1] * data_set[i][j]; //算bias以後的值
                    }
                    double denominator = 1 + Math.exp(-(decimal_weight[count_weight][0] * bias + num1));
                    hidden_ans[k] = 1 / denominator;
                }

                double num2 = 0;
                for(int k = 0 ; k < hidden_ans.length ; k++){
                    num2 += decimal_weight[count_weight][data_set[i].length*hidden_num + 1 + k] * hidden_ans[k];
                }

                double denominator = 1 + Math.exp(-(decimal_weight[count_weight][data_set[i].length*hidden_num] * bias + num2));
                double output_actual_ans = 1 / denominator;

                if((int) output_actual_ans*10 == (int)data_set[i][data_set[i].length -1]){
                    true_coount ++;
                }
            }

            fitness_original[count_weight] = Double.valueOf(String.format("%.02f", true_coount/(double)data_set.length));
        }
        return fitness_original;
    }

    //交配 function
    public static int [][] mating(int [][] binary_weight, int bytes_num){
        Random ran = new Random();
        int pointer_gene [] = new int [2];  //random 一個序列 1 ~ 10
        int [][] weight = new int[bytes_num][8];

        pointer_gene = randomPointer(2, 1, 10);

        for (int i = 0 ; i < pointer_gene.length ; i++){
            pointer_gene[i] = pointer_gene[i] - 1;
        }

        for(int i = 0 ; i < 5 ; i++){
            int num = (int) (0.1*(bytes_num*8-1)); // 有幾個1 (bits整數成上0.1)
            int amount_of_one = ran.nextInt(num)+1; //random
            int pointer_weight [] = new int [amount_of_one]; //random 位置

            pointer_weight = randomPointer(amount_of_one, 0, bytes_num*8-1);

            for(int j = 0 ; j < pointer_weight.length ; j++){
                int test = binary_weight[pointer_gene[0]][pointer_weight[j]];
                binary_weight[pointer_gene[0]][pointer_weight[j]] = binary_weight[pointer_gene[1]][pointer_weight[j]];
                binary_weight[pointer_gene[1]][pointer_weight[j]] = test;
            }
        }

        return binary_weight;
    }

    //突變 function
    public static int [][] mutation(int [][] binary_weight, int bytes_num){
        Random ran = new Random();

        for(int i = 0 ; i < binary_weight.length ; i++){
            int num = (int) (0.1*(bytes_num*8-1)); // 有幾個1 (bits整數成上0.1)
            int amount_of_one = ran.nextInt(num)+1; //random
            int pointer_weight [] = new int [amount_of_one]; //random 位置
            pointer_weight = randomPointer(amount_of_one, 0,bytes_num*8-1);

            for(int j = 0 ; j < pointer_weight.length ; j++){
                if(binary_weight[i][pointer_weight[j]] == 1)
                    binary_weight[i][pointer_weight[j]] = 0 ;
                if(binary_weight[i][pointer_weight[j]] == 0)
                    binary_weight[i][pointer_weight[j]] = 1 ;
            }
        }

        return binary_weight;
    }

    //挑選最好7條，與隨機3條
    public static  int [] bestOfFitness(double [] fitness_all){
        Random ran = new Random();

        int [] max_point = new int [10];
        for(int i = 0 ; i < max_point.length ; i++){
            max_point[i] = -99999;
        }

        for(int i = 0 ; i < 7 ; i++){
            double max = -99999;
            for(int j = 0 ; j < fitness_all.length ; j++){
                if(fitness_all[j] >= max){
                    int test = 0 ; //是否重覆
                    for(int k = 0 ; k < max_point.length ; k++){
                        if(max_point[k] == j)
                            test ++;
                    }
                    if(test == 0){
                        max = fitness_all[j];
                        max_point[i] = j;
                    }
                }
            }
        }

        for(int i = 7 ; i <= 9 ; i++){
            boolean uni = true;
            while (uni){
                int number = ran.nextInt(29)+0;
                int count = 0 ;

                for(int j = 0 ; j < max_point.length ; j++){
                    if(max_point[j] == number)
                        count++;
                }

                if(count == 0){
                    max_point[i] = number;
                    uni = false;
                }
            }
        }

        return max_point;
    }
}
