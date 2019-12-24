#include <iostream>
#include <string.h>
#include "fitsio.h"
#include <fstream>
#include <vector>
#include <math.h>

//prototype
void find_edges(std::vector<double>& x, std::vector<double>& y, std::vector<double> time_series);
void get_fake_elv(double RA, double DEC, std::vector<double> x, std::vector<double> y, std::vector<double> z, std::vector<double>& theta);
void PrintError(int status);
//
//
//


/* global variables */
int HV_criteria = -800;
double ra_source = 83.63322083;
double dec_source= 22.01446111111;
#define PI 3.14159265
double earch_theta = 70*PI/180;

int main(int argc, char* argv[])
{
    std::cout << "############# CsI Pulsation search #############" << std::endl;
    std::cout << ">>>>>>>>>>>>> Data reduction" << std::endl;
    char* filename = argv[1];
    char* hvfilename = argv[2];
    char* orbitname  = argv[3];
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
        if (highV[i] <= HV_criteria)
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
        for (int j=0; j<left_edges.size(); j++)
        {
            if (time_evt[i] >= left_edges[j] && time_evt[i] <= right_edges[j])
            {
                new_event.push_back(time_evt[i]);
                break;
            }
        }
    }
    /* Finish deleting GRB mode data */
    std::cout << std::fixed << "Entries of new Events = " << new_event.size() << std::endl;

    /* ############### calculate fake elv */
    // get column number
    int orbit_nRows, time_orbit_colnum, x_colnum, y_colnum, z_colnum;
    if (fits_open_file(&fptr, orbitname, READONLY, &status)) PrintError(status);
    if (ffmahd(fptr, 2, &datatype, &status)) PrintError(status);
    if (ffgky(fptr, TINT, "NAXIS2", &orbit_nRows, NULL, &status)) PrintError(status);
    if (fits_get_colnum(fptr, CASEINSEN, "TIME", &time_orbit_colnum, &status)) PrintError(status);
    if (fits_get_colnum(fptr, CASEINSEN, "X", &x_colnum, &status)) PrintError(status);
    if (fits_get_colnum(fptr, CASEINSEN, "Y", &y_colnum, &status)) PrintError(status);
    if (fits_get_colnum(fptr, CASEINSEN, "Z", &z_colnum, &status)) PrintError(status);
    std::cout << "nRows of Orbit = " << orbit_nRows << std::endl;

    double* tmp_time_orbit;
    double* tmp_j2000_x;
    double* tmp_j2000_y;
    double* tmp_j2000_z;
    tmp_time_orbit = new double[orbit_nRows];
    tmp_j2000_x = new double[orbit_nRows];
    tmp_j2000_y = new double[orbit_nRows];
    tmp_j2000_z = new double[orbit_nRows];
    if (fits_read_col(fptr, TDOUBLE, time_orbit_colnum, 1, 1, orbit_nRows, &doublenull, tmp_time_orbit, &anynull, &status)) PrintError(status);
    if (fits_read_col(fptr, TDOUBLE, x_colnum, 1, 1, orbit_nRows, &doublenull, tmp_j2000_x, &anynull, &status)) PrintError(status);
    if (fits_read_col(fptr, TDOUBLE, y_colnum, 1, 1, orbit_nRows, &doublenull, tmp_j2000_y, &anynull, &status)) PrintError(status);
    if (fits_read_col(fptr, TDOUBLE, z_colnum, 1, 1, orbit_nRows, &doublenull, tmp_j2000_z, &anynull, &status)) PrintError(status);
    if (fits_close_file(fptr, &status)) PrintError(status);
    std::vector<double> time_orbit;
    std::vector<double> j2000_x;
    std::vector<double> j2000_y;
    std::vector<double> j2000_z;
    for (int i=0; i<orbit_nRows; i++)
    {
        time_orbit.push_back(tmp_time_orbit[i]);
        j2000_x.push_back(tmp_j2000_x[i]);
        j2000_y.push_back(tmp_j2000_y[i]);
        j2000_z.push_back(tmp_j2000_z[i]);
    }

    // get the fake ELV angle
    std::vector<double> elv;
    std::vector<double> new_time_orbit;
    get_fake_elv(ra_source, dec_source, j2000_x, j2000_y, j2000_z, elv);
    for (int i=0; i<orbit_nRows; i++)
    {
        if (elv[i] >= earch_theta)
        {
            new_time_orbit.push_back(time_orbit[i]);
        }
    }

    if (new_time_orbit.size() != 0)
    {
        std::vector<double> left_edges_orbit;
        std::vector<double> right_edges_orbit;
        std::vector<double> new_event_orbit;
        find_edges(left_edges_orbit, right_edges_orbit, new_time_orbit);
        /* finish finding the GTI edges */
        /* select the event from the first GTI selection (new_event) by this good GTI */
        for (int i=0; i<new_event.size(); i++)
        {
            for (int j=0; j<left_edges.size(); j++)
            {
                if (new_event[i] >= left_edges_orbit[j] && new_event[i] <= right_edges_orbit[j])
                {
                    new_event_orbit.push_back(new_event[i]);
                    break;
                }
            }
        }
        std::cout << "THE FINAL ENTRIES = " << new_event_orbit.size() << " (" << new_event_orbit.size()*100/Evt_nRow << "%)" << std::endl;
    }
    else 
    {
        std::cout << "THE FINAL ENTRIES = " << 0 << std::endl;
    }


    // cleaning
    delete[] tmp_time_orbit, tmp_j2000_x, tmp_j2000_y, tmp_j2000_z;
    delete[] time_hv;
    delete[] highV;

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

void get_fake_elv(double RA, double DEC, std::vector<double> x, std::vector<double> y, std::vector<double> z, std::vector<double>& theta)
{

    double x_source = cos(DEC*PI/180) * cos(RA*PI/180);
    double y_source = cos(DEC*PI/180) * sin(RA*PI/180);
    double z_source = sin(DEC*PI/180);
    
    for (int i=0; i< x.size(); i++)
    {
        theta.push_back( acos((-x[i]*x_source - y[i]*y_source - z[i]*z_source)/(sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]))) );
        // the Theta value in units of Radius
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
