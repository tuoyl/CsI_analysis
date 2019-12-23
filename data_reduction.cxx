#include <iostream>
#include <string.h>
#include "fitsio.h"
#include <fstream>
#include <vector>

//prototype
void find_edges(std::vector<double>& x, std::vector<double>& y, std::vector<double> time_series);
void PrintError(int status);
//
//
//



int main(int argc, char* argv[])
{
    std::cout << "############# CsI Pulsation search #############" << std::endl;
    std::cout << ">>>>>>>>>>>>> Data reduction" << std::endl;
    char* filename = argv[1];
    char* hvfilename = argv[2];
    std::cout << filename << std::endl;
    std::cout << hvfilename << std::endl;

    /* Read the time and High Voltage column in HV file */
    fitsfile* fptr;
    int status = 0;
    int datatype, anynull;
    double doublenull = 0;
    int HV_nRows = 0;
    // get column number
    if (fits_open_file(&fptr, hvfilename, READONLY, &status)) PrintError(status);
    if (ffmahd(fptr, 2, &datatype, &status)) PrintError(status);
    if (ffgky(fptr, TINT, "NAXIS2", &HV_nRows, NULL, &status)) PrintError(status);
    std::cout << "nRows of HV = " << HV_nRows << std::endl;

    double* time_hv;
    double* highV;

    time_hv = new double[HV_nRows];
    highV   = new double[HV_nRows];

    int time_colnum;
    int hv_colnum;
    if (fits_get_colnum(fptr, CASEINSEN, "TIME", &time_colnum, &status)) PrintError(status);
    if (fits_get_colnum(fptr, CASEINSEN, "HV_PHODet_0", &hv_colnum, &status)) PrintError(status);
    if (fits_read_col(fptr, TDOUBLE, time_colnum, 1, 1, HV_nRows, &doublenull, time_hv, &anynull, &status)) PrintError(status);
    if (fits_read_col(fptr, TDOUBLE, hv_colnum, 1, 1, HV_nRows, &doublenull, highV, &anynull, &status)) PrintError(status);
    if (fits_close_file(fptr, &status)) PrintError(status);
    /* finish read HV file */


    /* Read the Event file */
    // get the column number of Time
    int Evt_nRow = 0;
    if (fits_open_file(&fptr, filename, READONLY, &status)) PrintError(status);
    if (ffmahd(fptr, 2, &datatype, &status)) PrintError(status);
    if (ffgky(fptr, TINT, "NAXIS2", &Evt_nRow, NULL, &status)) PrintError(status);
    std::cout << "nRows of Event = " << Evt_nRow << std::endl;

    double* time_evt;
    doublenull = 0;
    time_evt = new double[Evt_nRow];
    int evt_time_colnum;
    if (fits_get_colnum(fptr, CASEINSEN, "TIME", &evt_time_colnum, &status)) PrintError(status);
    if (fits_read_col(fptr, TDOUBLE, evt_time_colnum, 1, 1, Evt_nRow, &doublenull, time_evt, &anynull, &status)) PrintError(status);
    if (fits_close_file(fptr, &status)) PrintError(status);
    /* finish read Evt file */

    /* ################ DELETE grb mode time intervals */
    std::vector<double> new_event;

    /* find the GTI edges from HV time */
    std::vector<double> new_time_hv;

    for (int i=0; i<HV_nRows; i++)
    {
        if (highV[i] <= -800)
        {
            new_time_hv.push_back(time_hv[i]);
        }
    }
    std::vector<double> left_edges;
    std::vector<double> right_edges;
    find_edges(left_edges, right_edges, new_time_hv);
    /* finish finding the GTI edges */
    /* select the event in that good GTI */
    for (int i=0; i<Evt_nRow; i++)
    {
        for (int j=0; j<left_edges.size(); i++)
        {
            if (time_evt[i] >= left_edges[j] && time_evt[i] <= right_edges[j])
            {
                new_event.push_back(time_evt[i]);
                break;
            }
            continue;
        }
    }
    /* Finish deleting GRB mode data */
    std::cout << std::fixed << "Entries of new Events: " << new_event.size() << std::endl;

    /* ############### calculate fake elv */




    // cleaning
    delete time_hv;
    delete highV;

    return 0;
}


void find_edges(std::vector<double>& left_edges, std::vector<double>& right_edges, std::vector<double> time_series)
{
    /* find the continuous time edges for event series */

    std::vector<double> tmp_time_series;
    std::vector<double> tmp_time_series_roll;

    tmp_time_series = time_series;

    tmp_time_series_roll = time_series;

    // init the left edge
    left_edges.push_back(time_series.at(0));

    for (int i=0; i<tmp_time_series.size()-1; i++)
    {
        if (tmp_time_series_roll[i+1] - tmp_time_series[i] > 1)
        {
            right_edges.push_back(tmp_time_series[i]);
            left_edges.push_back(tmp_time_series_roll[i+1]);
        }
    }

    right_edges.push_back(time_series.back());
//

    for (int i=0; i<left_edges.size(); i++)
    {std::cout << "edges: " << std::fixed << left_edges[i] << " " << right_edges[i] << std::endl;
    }
}


void PrintError(int status)
{
    if (status)
    {
        fits_report_error(stderr, status);
        exit(status);
    }
    return;
}
