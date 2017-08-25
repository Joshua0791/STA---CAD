#include <map>
#include <algorithm>
#include <array>
#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <stdlib.h>
using namespace std;
class LUT {                   //Class for Look up table
public:
vector<string> Allgate_name; //all cells defined in the LUT
map<string,long double> cap; //map of cap values for a cell name of type string
map<string,long double**> num_in; //map of delays values for a cell name of type string
map<string,long double**> output_slew; //map of output slew values for a cell name of type string
void assignarrays(string); //function to pass the NLDM file name from which the above arrays can be populated
void alloc(int n) { Allgate_name.resize(n); } //resize vectors
map<string,long double*> tau_in; //corresponds to the 1st index in the LUT- check definition in the LUT template provided at the beginning of the LUT definition
map<string,long double*> cload; //corresponds to the 2nd index in the LUT
};

class node {                //Class for gates
    public:
        string name, outname; //name indicates cell type (nand, nor etc.), outname denotes the output wire name
        vector<string> input,out,in; //input names,primary outputs, primary inputs
        long double Cload; //load cap of this node
        vector<node*> inputs; //fanin nodes of this node
        vector<node*> outputs; //fanout nodes of this node
        vector<long double> Tau_in; //vector of input slews (for all inputs to the gate), to be used for STA
        vector<long double> inp_arrival; //vector of input arrival times for input transitions (ignore rise or fall)
        vector<long double> outp_arrival; //vector of output arrival times,outp_arrival= inp_arrival + cell_delay; cell_delay will be calculated from NLDM
        long double max_out_arrival; //arrival time at the output of this gate using max on (inp_arrival + cell_delay)
        long double Tau_out; //resulting output slew
        int num_in;         //number of inputs
        int processed_in;  // number of inputs processed
        int flag_out;       //flag to check of its a primary output
        int flag;  //processed or not
        long double rat,slack; // required arrival time, slack
        void outputs_alloc(int n) { outputs.resize(n); }    //resize vectors
        void alloc(int n) { input.resize(n); }  //resize vectors
        void inputs_alloc(int n) { inputs.resize(n); }  //resize vectors
};
LUT* lookup(ifstream &fp,ofstream &op1)
{
    LUT *p1;
    p1 = new LUT;
    int i=0,j=0;
    string s,stemp,stemp1,name,cellname;
    const char *pch ;
        if(fp.is_open()) {
                while(getline(fp, s)){
            if(s.find("cell ") < s.length()){
                p1->alloc(i+1);
                cellname=s.substr(s.find("(")+1, (s.find(")")-s.find("(")-1));
                stemp=s.substr(s.find("(")+1, (s.find(")")-s.find("(")-5));
                if(stemp.compare("IN")==0){stemp="NOT";}
                if(stemp.compare("BU")==0){stemp="BUFF";}
                p1->Allgate_name[i]=stemp;
                name=stemp;
                op1 << "Cell:" << cellname << "\n";
                i++;
                j=0;getline(fp,s);
                while(s.compare("")==0){getline(fp,s);
                }
                if(s.find("capacitance	") < s.length()){
                stemp1=(s.substr(s.find(":")+2 , 8));
                pch = stemp1.c_str();
                p1->cap[stemp]=(strtold(pch,NULL));
                op1 << "capacitance: " << p1->cap[stemp] << endl;
                getline(fp,s);
                getline(fp,s);
//index1
                j=0;
                stemp1=(s.substr(s.find("(")+2, (s.find(",") - s.find("(")-2)));
                pch = stemp1.c_str();
                p1->tau_in[name]= new long double [7];
                *(p1->tau_in[name])= strtold(pch,NULL);
                op1<< "Input slew: " << *(p1->tau_in[name]) << ",";
                stemp = (s.substr(s.find(",")+1,(s.find(")") - s.find(","))-1));
                j++;
                while(stemp.find(",") < stemp.length()-5) {
                    stemp1=(stemp.substr(0,stemp.find(",")));
                    pch = stemp1.c_str();
                    *(p1->tau_in[name]+j)= strtold(pch,NULL);
                    op1 << *(p1->tau_in[name]+j) << ",";
                    j++;
                    stemp = (stemp.substr(stemp.find(",")+1));
                }
                stemp1=(stemp.substr(0,(stemp.find(",")-1)));
                pch = stemp1.c_str();
                *(p1->tau_in[name]+j)= strtold(pch,NULL);
                op1 << *(p1->tau_in[name]+j) << ";" << endl;
                j=0;
//index2
                getline(fp,s);
                stemp1=(s.substr(s.find("(")+2, (s.find(",") - s.find("(")-2)));
                pch = stemp1.c_str();
                p1->cload[name]= new long double [7];
                *(p1->cload[name])= strtold(pch,NULL);
                op1<< "Load cap: " << *(p1->cload[name]) << ",";
                stemp = (s.substr(s.find(",")+1,(s.find(")") - s.find(","))-1));
                j++;
                while(stemp.find(",") < stemp.length()-5) {
                    stemp1=(stemp.substr(0,stemp.find(",")));
                    pch = stemp1.c_str();
                    *(p1->cload[name]+j)= strtold(pch,NULL);
                    op1 << *(p1->cload[name]+j) << ",";
                    j++;
                    stemp = (stemp.substr(stemp.find(",")+1));
                }
                stemp1=(stemp.substr(0,(stemp.find(",")-1)));
                pch = stemp1.c_str();
                *(p1->cload[name]+j)= strtold(pch,NULL);
                op1<<  *(p1->cload[name]+j) << ";" << endl;

// table
                getline(fp,s);
                if(s.find("values (") < s.length()){
                i=0;j=0;
                stemp1=(s.substr(s.find("(")+2, (s.find(",") - s.find("(")-2)));
                pch = stemp1.c_str();
                p1->num_in[name] = new long double*[7];
                for(int i = 0; i < 7; ++i)
                p1->num_in[name][i] = new long double[7];
                i=0;
                (p1->num_in[name])[i][j]= strtold(pch,NULL);
                op1<< "cell delays: \n" << *(p1->num_in[name]) << ",";
                stemp = (s.substr(s.find(",")+1,(s.find(")") - s.find(","))-1));
                j++;
                while(stemp.find(",") < stemp.length()-5) {
                    stemp1=(stemp.substr(0,stemp.find(",")));
                    pch = stemp1.c_str();
                    (p1->num_in[name])[i][j]= strtold(pch,NULL);
                    op1 << (p1->num_in[name])[i][j] << ",";
                    j++;
                    stemp = (stemp.substr(stemp.find(",")+1));
                    }
                stemp1=(stemp.substr(0,(stemp.find(",")-1)));
                pch = stemp1.c_str();
                p1->num_in[name][i][j]= strtold(pch,NULL);
                op1 << (p1->num_in[name])[i][j] << ";" << endl;
                j++;
                i++;
                while (i<=6){
                    j=0;
                    getline(fp,s);
                    stemp1=(s.substr(s.find("\"")+1, (s.find(",") - s.find("\"")-1)));
                    pch = stemp1.c_str();
                    p1->num_in[name][i][j]= strtold(pch,NULL);
                    op1 << p1->num_in[name][i][j] << ",";
                    stemp = (s.substr(s.find(",")+1,(s.find(")") - s.find(","))-1));
                    j++;
                    while(stemp.find(",") < stemp.length()-5) {
                        stemp1=(stemp.substr(0,stemp.find(",")));
                        pch = stemp1.c_str();
                        p1->num_in[name][i][j]= strtold(pch,NULL);
                        op1 << p1->num_in[name][i][j] << ",";
                        j++;
                        stemp = (stemp.substr(stemp.find(",")+1));
                        }
                    stemp1=(stemp.substr(0,(stemp.find(",")-1)));
                    pch = stemp1.c_str();
                    p1->num_in[name][i][j]= strtold(pch,NULL);
                    op1<< p1->num_in[name][i][j] << ";\n";
                    j++;
                    i++;
                    }
                getline(fp,s);
                getline(fp,s);
                while(s.compare("")==0){getline(fp,s);
                }
                getline(fp,s);
                getline(fp,s);
                getline(fp,s);
                if(s.find("values (") < s.length()){
                i=0;j=0;
                p1->output_slew[name] = new long double*[7];
                for(int i = 0; i < 7; ++i)
                p1->output_slew[name][i] = new long double[7];
                i=0;
                stemp1=(s.substr(s.find("(")+2, (s.find(",") - s.find("(")-2)));
                pch = stemp1.c_str();
                p1->output_slew[name][i][j]= strtold(pch,NULL);
                                op1<< "output slews: \n" << p1->output_slew[name][i][j] << ",";
                stemp = (s.substr(s.find(",")+1,(s.find(")") - s.find(","))-1));
                                j++;
                while(stemp.find(",") < stemp.length()-5) {
                                    stemp1=(stemp.substr(0,stemp.find(",")));
                    pch = stemp1.c_str();
                    p1->output_slew[name][i][j]= strtold(pch,NULL);
                                    op1 << p1->output_slew[name][i][j] << ",";
                                    j++;
                    stemp = (stemp.substr(stemp.find(",")+1));
                    }
                stemp1=(stemp.substr(0,(stemp.find(",")-1)));
                pch = stemp1.c_str();
                p1->output_slew[name][i][j]= strtold(pch,NULL);
                                op1 << p1->output_slew[name][i][j] << ";" << endl;
                        j++;
                i++;
                while (i<=6){
                    j=0;
                    getline(fp,s);
                    stemp1=(s.substr(s.find("\"")+1, (s.find(",") - s.find("\"")-1)));
                    pch = stemp1.c_str();
                    p1->output_slew[name][i][j]= strtold(pch,NULL);
                                    op1 << p1->output_slew[name][i][j] << ",";
                    stemp = (s.substr(s.find(",")+1,(s.find(")") - s.find(","))-1));
                                j++;
                    while(stemp.find(",") < stemp.length()-5) {
                                        stemp1=(stemp.substr(0,stemp.find(",")));
                        pch = stemp1.c_str();
                        p1->output_slew[name][i][j]= strtold(pch,NULL);
                                        op1 << p1->output_slew[name][i][j] << ",";
                                        j++;
                        stemp = (stemp.substr(stemp.find(",")+1));
                        }
                    stemp1=(stemp.substr(0,(stemp.find(",")-1)));
                    pch = stemp1.c_str();
                    p1->output_slew[name][i][j]= strtold(pch,NULL);
                                    op1<< p1->output_slew[name][i][j] << ";\n";
                    j++;
                    i++;
                    }
                    op1 << "------------------------------------------------------------------------" <<endl;


                }

                }
                }

            }
        }
    }
    op1.close();
    return p1;
}

vector<node*> netlist(ifstream &fp,ofstream &op2)
{
    vector<node*> gates;
    vector<node*> que;
    vector<string> out,in;
    node *p1;
    int gate_cnt=0,out_cnt=0,in_cnt=0,nand_cnt=0,nor_cnt=0,or_cnt=0,and_cnt=0,not_cnt=0,xor_cnt=0,buff_cnt=0;
    string s,stemp;
    if(fp.is_open()) {
        while(getline(fp, s)){
            if(s.find("=") < s.length()){
                p1= new node ;
                (*p1).num_in=0;
                p1->flag_out=0;
                p1->outname=s.substr(0, s.find("=")-1);
                (*p1).name=s.substr(s.find("=")+2,(s.find("(")-s.find("=")-2));
                stemp=p1->name;
                if(stemp.compare("NAND")==0){ nand_cnt++;}
                if(stemp.compare("NOR")==0){ nor_cnt++;}
                if(stemp.compare("OR")==0){ or_cnt++;}
                if(stemp.compare("AND")==0){ and_cnt++;}
                if(stemp.compare("BUFF")==0) { buff_cnt++;}
                if(stemp.compare("XOR")==0) { xor_cnt++;}
                if(stemp.compare("NOT")==0) { not_cnt++;}
                if(stemp.compare("NOT")==0 || stemp.compare("BUFF")==0) {
                            p1 -> alloc((*p1).num_in+1);
                            (*p1).input[(*p1).num_in]=(s.substr(s.find("(")+1, (s.find(")") - s.find("("))-1));
                            (*p1).num_in++;
                            gate_cnt++;
                            gates.resize(gate_cnt);
                            gates[gate_cnt-1]=p1;
                }
                else {
                p1 -> alloc((*p1).num_in+1);
                (*p1).input[(*p1).num_in]=(s.substr(s.find("(")+1, (s.find(",") - s.find("("))-1));
                stemp = s.substr(s.find(",")+2,(s.find(")") - s.find(","))-2);
                (*p1).num_in++;
                while(stemp.find(",") < stemp.length()) {
                    p1 -> alloc((*p1).num_in+1);
                    (*p1).input[(*p1).num_in]=(stemp.substr(0,stemp.find(",")));
                    (*p1).num_in++;
                    stemp = stemp.substr(stemp.find(",")+2);
                }
                p1 -> alloc((*p1).num_in+1);
                (*p1).input[(*p1).num_in]=(stemp);
                (*p1).num_in++;
                gate_cnt++;

                gates.resize(gate_cnt);
                gates[gate_cnt-1]=p1;
                }
            }
            if(s.find("OUTPUT") < s.length()){
                out.resize(out_cnt+1);
                out[out_cnt]=s.substr(s.find("(")+1, (s.find(")") - s.find("("))-1);
                out_cnt++;
            }
            if(s.find("INPUT") < s.length()){
                in.resize(in_cnt+1);
                in[in_cnt]=s.substr(s.find("(")+1, (s.find(")") - s.find("("))-1);
                in_cnt++;
            }

        }
}
    op2 << in_cnt << " primary inputs\n";
    op2 << out_cnt << " primary outputs\n";
    op2 << nand_cnt << " NAND gates\n";
    op2 << nor_cnt << " NOR gates\n";
    op2 << or_cnt << " OR gates\n";
    op2 << and_cnt << " AND gates\n";
    op2 << not_cnt << " NOT gates\n";
    op2 << buff_cnt << " BUFF gates\n";
    op2 << xor_cnt << " XOR gates\n";
    op2 << gate_cnt << " Total gates\n";
    op2 << "\nFANOUT \n";
    int i=0,j=0,l=0,k=0;

    while(i<gates.size()){
        l=0;
        j=0;
        node *a1=gates[i];
        string s1=a1->outname;
        op2 << a1->name << "-" << s1 <<":";
        while(j<gates.size()){
            node *a2=gates[j];
            for(k=0;k<(*a2).num_in;k++){
                string s2=a2->input[k];
                if((s1.compare(s2))==0) {
                    a1->outputs.push_back(a2);
                    op2<< " " << a2->name<<"-" << (*a2).outname ;
                    l++;
                }
            }
            j++;
        }
        j=0;
        while(j<out_cnt){
            if((s1.compare(out[j])==0)) {
                a1->flag_out++;
                a1->out.push_back(out[j]);
                op2 << " OUTP";
            }
            j++;
        }
        i++;
        op2 << "\n" ;
    }
    op2 << "\nFANIN \n";
    i=0;j=0;l=0;k=0;
    int flag=0;
    while(i<gates.size()){
        l=0;
        j=0;
        string s1;
        node *a1=gates[i];
        a1->flag=0;
        a1->processed_in=0;
        a1->rat=0;
        flag=0;
        op2 << a1->name << "-" << a1->outname <<":" ;
        for(k=0;k<(*a1).num_in;k++){
            s1=a1->input[k];
            j=0;
            while(j<gates.size()){
                node *a2=gates[j];
                string s2=a2->outname;
                if((s1.compare(s2))==0) {
                    a1->inputs.push_back(a2);
                    op2<< " " << a2->name<<"-" << a2->outname << "" ;
                    l++;
                }
                j++;
            }
            j=0;
            while(j<in_cnt){
                if((s1.compare(in[j])==0)) {
                    flag++;
                    a1->processed_in++;

                    a1->in.push_back(in[j]);
                    if( flag == a1->num_in ) {
                            que.push_back(a1);}
                    op2 << " INP-" << s1;
                }
                j++;
            }
        }
        i++;
        op2 << "\n";
    }
op2.close();
return que;
}

long double find_delay (long double tau_in, long double loadcap, string name, int fanin, LUT* m1)
{
   long double final_delay=0;
   int i=0,j=0,flag=-1,flag1=-1,t1,t2,c1,c2;
   long double slew1,slew2,cap2,cap1;
    for(i=0;i<6;i++){
        if(tau_in < m1->tau_in[name][i]){ flag=i; break;}
    }
    if(flag<=6 && flag>0){t1=flag-1;t2=flag;slew1=m1->tau_in[name][flag-1];slew2=m1->tau_in[name][flag];}
    if(flag==0){ t1=0;t2=1;slew1=m1->tau_in[name][0];slew2=m1->tau_in[name][1];}
    if(flag==-1){t1=5;t2=6;slew1=m1->tau_in[name][5];slew2=m1->tau_in[name][6];}
    for(i=0;i<6;i++){
        if(loadcap < m1->cload[name][i]){ flag1=i; break;}
    }
    if(flag1<=6 && flag1>0){c1=flag1-1;c2=flag1;cap1=m1->cload[name][flag1-1];cap2=m1->cload[name][flag1];}
    if(flag1==0){ c1=0;c2=1;cap1=m1->cload[name][0];cap2=m1->cload[name][1];}

    if(flag1==-1){c1=5;c2=6;cap1=m1->cload[name][5];cap2=m1->cload[name][6];}
    final_delay=(m1->num_in[name][t1][c1])*(cap2-loadcap)*(slew2-tau_in);
    final_delay+=(m1->num_in[name][t1][c2])*(loadcap-cap1)*(slew2-tau_in);
    final_delay+=(m1->num_in[name][t2][c1])*(cap2-loadcap)*(tau_in-slew1);
    final_delay+=(m1->num_in[name][t2][c2])*(loadcap-cap1)*(tau_in-slew1);
    final_delay=final_delay/((cap2- cap1)*(slew2-slew1));
    if(fanin>2){ final_delay=(fanin)/2.0*final_delay;  }
    if(final_delay<0) { final_delay=0;}
    return final_delay;

}

long double find_slew (long double tau_in, long double loadcap, string name, int fanin, LUT* m1)
{
   long double final_slew=0;
   int i=0,j=0,flag=-1,flag1=-1,t1,t2,c1,c2;
   long double slew1,slew2,cap2,cap1;
    for(i=0;i<6;i++){
        if(tau_in < m1->tau_in[name][i]){ flag=i; break;}
    }
    if(flag<=6 && flag>0){t1=flag-1;t2=flag;slew1=m1->tau_in[name][flag-1];slew2=m1->tau_in[name][flag];}
    if(flag==0){ t1=0;t2=1;slew1=m1->tau_in[name][0];slew2=m1->tau_in[name][1];}
    if(flag==-1){t1=5;t2=6;slew1=m1->tau_in[name][5];slew2=m1->tau_in[name][6];}
    for(i=0;i<6;i++){
        if(loadcap < m1->cload[name][i]){ flag1=i; break;}}
    if(flag1<=6 && flag1>0){c1=flag1-1;c2=flag1;cap1=m1->cload[name][flag1-1];cap2=m1->cload[name][flag1];}
    if(flag1==0){ c1=0;c2=1;cap1=m1->cload[name][0];cap2=m1->cload[name][1];}
    if(flag1==-1){c1=5;c2=6;cap1=m1->cload[name][5];cap2=m1->cload[name][6];}
    final_slew=(m1->output_slew[name][t1][c1])*(cap2-loadcap)*(slew2-tau_in);
    final_slew+=(m1->output_slew[name][t1][c2])*(loadcap-cap1)*(slew2-tau_in);
    final_slew+=(m1->output_slew[name][t2][c1])*(cap2-loadcap)*(tau_in-slew1);
    final_slew+=(m1->output_slew[name][t2][c2])*(loadcap-cap1)*(tau_in-slew1);
    final_slew=final_slew/((cap2- cap1)*(slew2-slew1));
    if(fanin>2){ final_slew=(fanin)/2.0*final_slew;  }
    if(final_slew<0) { final_slew=0.002;}
    return final_slew;

}


int main(int argc, char *argv[])
{
    string path1,path2,opath1,opath2;
    vector<string> prim_in,critical_path;
    vector<node*> que;
    vector<node*> que1;
    que1.push_back(0);
    map<string,long double> inp_rat,inp_cp,out_rat,out_tau;
    node* p1;
    node* p2;
    node* p3;
    LUT* m1;
    int j=0,i=0,k=0,l=0,m=0,flag=0;
    long double max_arrival=0,temp=0,min_slack=0;
    path1=argv[1];
    ifstream InputFile1 (path1.c_str());
    path2=argv[2];
    ifstream InputFile2 (path2.c_str());
    opath1=argv[3];
    opath2=argv[4];
    ofstream OutputFile1;
    OutputFile1.open(opath1.c_str());
    ofstream OutputFile2;
    OutputFile2.open(opath2.c_str());
    m1=lookup(InputFile1,OutputFile1);
    que=netlist(InputFile2,OutputFile2);
    while( j < que.size()  ){
        i++;
        p1=que[j];
        k=0;
        p1->Cload=0;
        if (p1->flag_out==1)
        {
          p1->Cload+=4*1.700230;
        }
        while (k < p1->outputs.size()) {
            p2 = p1->outputs[k];
            p1->Cload+=m1->cap[p2->name];
            p2->processed_in++;
            if(p2->processed_in == p2->num_in && p2->flag == 0){
                p2->flag=1;
                que.push_back(p2);
            }
            k++;
        }
        k=0;
        long double delay;
        if(p1->inputs.size() < p1->num_in){
            for(m=0;m<(p1->num_in - p1->inputs.size());m++){
            p1->Tau_in.push_back(0.002);
            p1->inp_arrival.push_back(0);
            delay= find_delay(p1->Tau_in[k],p1->Cload,p1->name,p1->num_in,m1);
            temp=(p1->inp_arrival[k]+delay);
            p1->outp_arrival.push_back(temp);
            k++;
            }
        }
        l=0;
        while (l < p1->inputs.size()) {
            p3 = p1->inputs[l];
            p1->Tau_in.push_back(p3->Tau_out);
            p1->inp_arrival.push_back(p3->max_out_arrival);
            delay= find_delay(p1->Tau_in[k],p1->Cload,p1->name,p1->num_in,m1);
            p1->outp_arrival.push_back((p1->inp_arrival[k]+delay));
            k++;
            l++;
        }

        p1->max_out_arrival=0;
        for(l=0;l<p1->outp_arrival.size();l++){
              if (p1->max_out_arrival < p1->outp_arrival[l]){
                p1->max_out_arrival= p1->outp_arrival[l];
                k=l;
              }

        }
        p1->Tau_out=find_slew(p1->Tau_in[k],p1->Cload,p1->name,p1->num_in,m1);
        if(p1->max_out_arrival>max_arrival && p1->flag_out==1)
        {
            max_arrival= p1->max_out_arrival;
            que1[0]=p1;
        }
        j++;

    }
cout << endl << "########################################################################################################################" << endl;
cout << endl << "Max delay=" << max_arrival*1000 << "ps" << endl;
long double rat_final=1.1*max_arrival;
min_slack=max_arrival;
i=0;j=0;k=0;l=0;
j=que.size()-1;
while( i < que.size()  ){
    p1=que[j];
    if(p1->out.size()!=0){
        out_rat[p1->out[0]]=p1->max_out_arrival;
        out_tau[p1->out[0]]=p1->Tau_out;
        if(p1->rat==0){p1->rat=rat_final; p1->slack=rat_final-p1->max_out_arrival;
        }
        else {if(p1->rat > rat_final){p1->rat=rat_final;p1->slack=rat_final-p1->max_out_arrival;}    }}
    for(k=0;k<p1->inputs.size();k++) {
       p2 = p1->inputs[k];
       if(p2->rat==0 || p2->rat > (p1->rat-find_delay(p1->Tau_in[k+p1->in.size()],p1->Cload,p1->name,p1->num_in,m1))){
        p2->rat=p1->rat-find_delay(p1->Tau_in[k+p1->in.size()],p1->Cload,p1->name,p1->num_in,m1);
        p2->slack=p2->rat -p2->max_out_arrival;
       }
    }
    l=0;
    while(p1->in.size()>l)
    {
        if(inp_rat[p1->in[l]]==0){inp_rat[p1->in[l]]=p1->rat-(p1->outp_arrival[l]-p1->inp_arrival[l]);}
        else {if( inp_rat[p1->in[l]] > (p1->rat-(p1->outp_arrival[l]-p1->inp_arrival[l]))) {
        inp_rat[p1->in[l]]=p1->rat-(p1->outp_arrival[l]-p1->inp_arrival[l]);}}
        inp_cp[p1->in[l]]+=m1->cap[p1->name];
        l++;
    }

    j--;
    i++;
}
k=0;
critical_path.push_back("OUTPUT-" + que1[0]->outname);
critical_path.push_back(que1[0]->name + "-" + que1[0]->outname);
while( k < que1.size()  )
{
   critical_path.push_back("0");
   p1=que1[k];
   if(p1->in.size() > 0) { p3=p1;}
   if(p1->inputs.size() > 0) {que1.push_back(p1->inputs[0]); critical_path[k+2]=(p1->inputs[0]->name + "-" + p1->inputs[0]->outname);
   path1=("Input-" + p1->inputs[0]->outname);
   if(p1->in.size() > 0) { if(p1->outname.compare("676")==0){cout << "here";}p3=p1;}
   min_slack=p1->inputs[0]->slack;
   }
   for(j=0;j<p1->inputs.size();j++)
   {
       p2=p1->inputs[j];
       if(p2->slack<min_slack)
       {
          min_slack=p2->slack;
          critical_path[k+2]= (p2->name + "-" + p2->outname);
          path1=("Input-" + p2->outname);
          if(p2->in.size() > 0) {p3=p2;}
          que1[k+1]=p2;
       }

   }
   k++;
}
cout << endl << "########################################################################################################################" << endl;
cout << endl << "Gate\t\tLoadCap(fF}\tArrival times(ps)\tSlews(ps)\tSlack(ps)\n";
map<string,long double>::const_iterator it = inp_cp.cbegin();
for(map<string,long double>::const_iterator klem=inp_rat.begin() ; klem!=inp_rat.cend() ;klem++)
{
   std::cout << "INPUT-" <<klem->first << "\t" << it->second << "\t\t0\t\t\t2\t\t" <<  klem->second*1000 << "\n";
   it++;
}
it = out_tau.cbegin();
for(map<string,long double>::const_iterator klem=out_rat.begin() ; klem!=out_rat.cend() ;klem++)
{
   std::cout << "OUTPUT("<<klem->first << ")\t0\t\t" << klem->second*1000 << "\t\t\t" << it->second*1000 << "\t\t" << klem->second*100 << "\n";
}
i=0;
while( i < que.size()  ){
    p1=que[i];
    cout<< p1->name<< "-" << p1->outname <<"\t"<< p1->Cload <<"\t\t"<< p1->max_out_arrival*1000 <<"\t\t\t"<< p1->Tau_out*1000 << "\t\t"<<p1->slack*1000 << "\n";
    i++;

}

i=0;
cout << endl << "########################################################################################################################" << endl;
cout << "Critical Path:" << endl;
if(p3->in.size() > 0){cout << "Input-" << p3->in[0];}
while( i < critical_path.size()-1){
    cout  << "->" <<critical_path[critical_path.size()-2-i] ;
    i++;
}
cout << endl <<"########################################################################################################################" << endl;

}


