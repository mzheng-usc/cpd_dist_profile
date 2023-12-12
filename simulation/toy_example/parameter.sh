#!/bin/bash
for a in `seq 0 0.1 1`; 
  do 
    for b in {1..50};
      do
        echo "$a $b"
      done;
  done;

#type_seq<-seq(0,1,0.1)
#type_seq<-seq(0,0.4,0.04)
#type_seq<-seq(2,22,2)
#type_seq<-seq(0.0000000001,0.5000000001,0.05)
