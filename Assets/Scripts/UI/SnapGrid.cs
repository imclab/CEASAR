﻿using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.SceneManagement;

public class SnapGrid : MonoBehaviour
{
    public GameObject snapItemPrefab;
    private List<GameObject> snaps;

    private SimulationManager manager { get {return SimulationManager.GetInstance();}}

    private void Start()
    {
        snaps = new List<GameObject>();
    }
    
    public void AddSnapItem(Pushpin newSnap)
    {
        if (snaps == null) snaps = new List<GameObject>();
        string snapText = newSnap.LocationName + ":\n" + newSnap.SelectedDateTime.ToShortDateString() + " " + newSnap.SelectedDateTime.ToShortTimeString();
        GameObject snapItem = (GameObject)Instantiate(snapItemPrefab, transform);
        snapItem.GetComponent<SnapItem>().snapItemText.text = snapText;
        snapItem.GetComponent<SnapItem>().snapshot = newSnap;
        snaps.Add(snapItem);
    }
}
