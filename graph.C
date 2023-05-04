void graph(){
    TCanvas* c = new TCanvas("c", "c", 800, 600);
    TMultiGraph* multiG = new TMultiGraph();
    TGraph* n = new TGraph("puntos.txt");
    TF1* a = new TF1("num", "1.1392070132*x - 0.03920701320/(x*x) - 0.3*sin(TMath::Log(x)) - 0.1*cos(TMath::Log(x))", 1,2);
    TGraph* aa = new TGraph(a);
    n->SetLineWidth(5);
    a->SetLineWidth(5);

    multiG->Add(n, "AC");
    multiG->Add(aa, "AC");
    multiG->Draw("A");

    TLegend* l = new TLegend(.1,.6,.5,.9,"");
    l->AddEntry(n, "Sol. Numerica", "LP");
    l->AddEntry(aa, "Sol. Analitica", "LP");
    l->Draw();


    c->Print("graph.png");
}