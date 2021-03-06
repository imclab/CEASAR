﻿using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using UnityEngine;

public enum LogLevel { Verbose, Info, Warning, Error }
public enum LogMessageCategory {All, Rendering, Networking, VR, UI, Interaction, Event, EventLog}

public static class CCDebug
{
    public static LogLevel CurrentLevel = LogLevel.Info;

    public static LogMessageCategory[] Categories = new[]
    {
        LogMessageCategory.All,
        LogMessageCategory.Event,
        LogMessageCategory.Interaction,
        LogMessageCategory.Networking,
        LogMessageCategory.Rendering,
        LogMessageCategory.UI,
        LogMessageCategory.VR,
        LogMessageCategory.EventLog
    };
    
    public static void Log(object logMessage, LogLevel level, LogMessageCategory category)
    {
        if (CurrentLevel <= level && Categories.Contains(category))
        {
            switch (level)
            {
                case LogLevel.Info:
                    Debug.Log(logMessage);
                    break;
                case LogLevel.Verbose:
                    Debug.Log(logMessage);
                    break; 
                case LogLevel.Warning:
                    Debug.LogWarning(logMessage);
                    break;
                case LogLevel.Error:
                    Debug.LogError(logMessage);
                    break;
            }
        }
    }

    public static void Log(object logMessage)
    {
        Log(logMessage, LogLevel.Verbose, LogMessageCategory.All);
    }
    public static void LogError(object logMessage)
    {
        Log(logMessage, LogLevel.Error, LogMessageCategory.All);
    }
    public static void LogWarning(object logMessage)
    {
        Log(logMessage, LogLevel.Warning, LogMessageCategory.All);
    }
    
}
