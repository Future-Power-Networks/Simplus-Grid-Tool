clear all
clc
close all

ListBus = [1,1;
           2,1];
       

       t1 = [4,4;
            2,2;
            NaN,0;
            1,1;
            NaN,0];
        
       t2 =  sortrows(t1,1)