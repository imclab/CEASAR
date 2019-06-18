﻿using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class ConstellationsController : MonoBehaviour
{
    List<Constellation> constellations = new List<Constellation>();
    public float lineWidth = .035f;

    public void AddConstellation(Constellation constellation)
    {
        constellations.Add(constellation);
    }

    public Constellation GetConstellation(string cFullName)
    {
        int index = constellations.FindIndex(el => el.constellationNameAbbr == cFullName);
        if (index < 0)
        {
            return constellations[index];
        }
        else
        {
            return null;
        }
    }

    public void HighlightSingleConstellation(string cFullName)
    {
        foreach (Constellation constellation in constellations)
        {
            constellation.Highlight(constellation.constellationNameFull == cFullName);
            constellation.ShowConstellationLines(constellation.constellationNameFull == cFullName, lineWidth);
        }
    }

    public void HighlightAllConstellations(bool highlight)
    {
        foreach (Constellation constellation in constellations)
        {
            constellation.Highlight(highlight);
            constellation.ShowConstellationLines(highlight, lineWidth);
        }
    }

    public void ShowSingleConstellation(string cFullName)
    {
        foreach (Constellation constellation in constellations)
        {
            constellation.ShowConstellationLines(constellation.constellationNameFull == cFullName, lineWidth);
        }
    }

    public void ShowAllConstellations(bool show)
    {
        foreach (Constellation constellation in constellations)
        {
            constellation.ShowConstellationLines(show, lineWidth);
        }
    }

}