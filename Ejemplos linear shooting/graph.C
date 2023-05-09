void graph(){
    TCanvas* c = new TCanvas("c", "c", 800, 600);
    TMultiGraph* multiG = new TMultiGraph();
    TGraph* n = new TGraph("puntos.txt");
    TF1* a = new TF1("num", "-(sin(x) + 3*cos(x))/10", 0,3.14159265/2);
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

    multiG->SetTitle("y'' = y' + 2y + cos(x)");
    c->Print("Exe3.png");
}