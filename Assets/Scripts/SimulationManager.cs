﻿using System.Collections.Generic;
using System.Text;
using UnityEngine;

// This class is the manager Singleton, and contains specific references to application-level objects
public class SimulationManager
{

    // can't use constructor, guaranteed singleton
    protected SimulationManager() { }
    private static SimulationManager instance;

    public static SimulationManager GetInstance()
    {
        return instance ?? (instance = new SimulationManager());
    }
    public System.Random rng = new System.Random();

    public GameObject NetworkControllerObject;

    private GameObject celestialSphereObject;

    public GameObject CelestialSphereObject
    {
        get { return celestialSphereObject; }
        set
        {
            celestialSphereObject = value;
            DataControllerComponent = celestialSphereObject.GetComponent<DataController>();
            MarkersControllerComponent = celestialSphereObject.GetComponentInChildren<MarkersController>();
            ConstellationsControllerComponent = celestialSphereObject.GetComponentInChildren<ConstellationsController>();
        }
    }

    public DataController DataControllerComponent { get; private set; }
    public MarkersController MarkersControllerComponent { get; private set; }
    public ConstellationsController ConstellationsControllerComponent { get; private set; }

    public bool IsReady = false;

    public string[] AnimalNames;
    public List<string> ColorNames = new List<string>();
    public List<Color> ColorValues = new List<Color>();

    public Color LocalPlayerColor = Color.white;

    // initial setup scale
    public readonly float InitialRadius = 100;
    public float CurrentScaleFactor(float sceneRadius)
    {
        return sceneRadius / InitialRadius;
    }

    // Random color (capitalized), random animal (capitalized), random number
    public string GenerateUsername()
    {
        int colorIndex = rng.Next(ColorNames.Count - 1);
        int animalIndex = rng.Next(AnimalNames.Length - 1);
        string randomNumber = rng.Next(999).ToString();
        return ColorNames[colorIndex].FirstCharToUpper() + AnimalNames[animalIndex].FirstCharToUpper() + randomNumber;
    }

    // We can find out the color value from the username
    public Color GetColorForUsername(string name)
    {
        StringBuilder sb = new StringBuilder();
        bool found = false;
        int i = 0;
        while (!found && i < name.Length)
        {
            if (char.IsUpper(name[i]) && i > 1)
            {
                found = true;
            }
            else
            {
                sb.Append(name[i]);
                i++;
            }
        }

        string colorName = sb.ToString();
        // Don't forget to lowercase the name!
        if (ColorNames.Contains(colorName.ToLower()))
        {
            return ColorValues[ColorNames.IndexOf(colorName.ToLower())];
        }
        else
        {
            Debug.Log("Color not found for " + colorName.ToLower() + " as part of " + name);
            return UnityEngine.Random.ColorHSV(0f, 1f, 1f, 1f, 0.9f, 1f);
        }
    }

    public float GetRelativeMagnitude(float starMagnitude)
    {
        //float min = DataControllerComponent.minMag;
        float max = DataControllerComponent.maxMag + Mathf.Abs(DataControllerComponent.minMag);
        return max - (starMagnitude + Mathf.Abs(DataControllerComponent.minMag));
    }
}