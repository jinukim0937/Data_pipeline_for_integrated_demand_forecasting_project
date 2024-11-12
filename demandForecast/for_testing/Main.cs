/*
* MATLAB Compiler: 8.4 (R2022a)
* Date: Fri Mar  3 15:29:16 2023
* Arguments:
* "-B""macro_default""-W""dotnet:demandForecast,Main,4.0,private,version=1.0""-T""link:lib
* ""-d""C:\Users\TesT\OneDrive - 연세대학교 (Yonsei University)\바탕
* 화면\YG1\demandForecast\for_testing""-v""class{Main:C:\Users\TesT\OneDrive -
* 연세대학교 (Yonsei University)\바탕 화면\YG1\demandForecast.m}"
*/
using System;
using System.Reflection;
using System.IO;
using MathWorks.MATLAB.NET.Arrays;
using MathWorks.MATLAB.NET.Utility;

#if SHARED
[assembly: System.Reflection.AssemblyKeyFile(@"")]
#endif

namespace demandForecast
{

  /// <summary>
  /// The Main class provides a CLS compliant, MWArray interface to the MATLAB functions
  /// contained in the files:
  /// <newpara></newpara>
  /// C:\Users\TesT\OneDrive - 연세대학교 (Yonsei University)\바탕
  /// 화면\YG1\demandForecast.m
  /// </summary>
  /// <remarks>
  /// @Version 1.0
  /// </remarks>
  public class Main : IDisposable
  {
    #region Constructors

    /// <summary internal= "true">
    /// The static constructor instantiates and initializes the MATLAB Runtime instance.
    /// </summary>
    static Main()
    {
      if (MWMCR.MCRAppInitialized)
      {
        try
        {
          Assembly assembly= Assembly.GetExecutingAssembly();

          string ctfFilePath= assembly.Location;

		  int lastDelimiter = ctfFilePath.LastIndexOf(@"/");

	      if (lastDelimiter == -1)
		  {
		    lastDelimiter = ctfFilePath.LastIndexOf(@"\");
		  }

          ctfFilePath= ctfFilePath.Remove(lastDelimiter, (ctfFilePath.Length - lastDelimiter));

          string ctfFileName = "demandForecast.ctf";

          Stream embeddedCtfStream = null;

          String[] resourceStrings = assembly.GetManifestResourceNames();

          foreach (String name in resourceStrings)
          {
            if (name.Contains(ctfFileName))
            {
              embeddedCtfStream = assembly.GetManifestResourceStream(name);
              break;
            }
          }
          mcr= new MWMCR("",
                         ctfFilePath, embeddedCtfStream, true);
        }
        catch(Exception ex)
        {
          ex_ = new Exception("MWArray assembly failed to be initialized", ex);
        }
      }
      else
      {
        ex_ = new ApplicationException("MWArray assembly could not be initialized");
      }
    }


    /// <summary>
    /// Constructs a new instance of the Main class.
    /// </summary>
    public Main()
    {
      if(ex_ != null)
      {
        throw ex_;
      }
    }


    #endregion Constructors

    #region Finalize

    /// <summary internal= "true">
    /// Class destructor called by the CLR garbage collector.
    /// </summary>
    ~Main()
    {
      Dispose(false);
    }


    /// <summary>
    /// Frees the native resources associated with this object
    /// </summary>
    public void Dispose()
    {
      Dispose(true);

      GC.SuppressFinalize(this);
    }


    /// <summary internal= "true">
    /// Internal dispose function
    /// </summary>
    protected virtual void Dispose(bool disposing)
    {
      if (!disposed)
      {
        disposed= true;

        if (disposing)
        {
          // Free managed resources;
        }

        // Free native resources
      }
    }


    #endregion Finalize

    #region Methods

    /// <summary>
    /// Provides a single output, 0-input MWArrayinterface to the demandForecast MATLAB
    /// function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// parseObj.addParameter('numForecasts',3);
    /// </remarks>
    /// <returns>An MWArray containing the first output argument.</returns>
    ///
    public MWArray demandForecast()
    {
      return mcr.EvaluateFunction("demandForecast", new MWArray[]{});
    }


    /// <summary>
    /// Provides a single output, 1-input MWArrayinterface to the demandForecast MATLAB
    /// function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// parseObj.addParameter('numForecasts',3);
    /// </remarks>
    /// <param name="dbName">Input argument #1</param>
    /// <returns>An MWArray containing the first output argument.</returns>
    ///
    public MWArray demandForecast(MWArray dbName)
    {
      return mcr.EvaluateFunction("demandForecast", dbName);
    }


    /// <summary>
    /// Provides a single output, 2-input MWArrayinterface to the demandForecast MATLAB
    /// function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// parseObj.addParameter('numForecasts',3);
    /// </remarks>
    /// <param name="dbName">Input argument #1</param>
    /// <param name="dbID">Input argument #2</param>
    /// <returns>An MWArray containing the first output argument.</returns>
    ///
    public MWArray demandForecast(MWArray dbName, MWArray dbID)
    {
      return mcr.EvaluateFunction("demandForecast", dbName, dbID);
    }


    /// <summary>
    /// Provides a single output, 3-input MWArrayinterface to the demandForecast MATLAB
    /// function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// parseObj.addParameter('numForecasts',3);
    /// </remarks>
    /// <param name="dbName">Input argument #1</param>
    /// <param name="dbID">Input argument #2</param>
    /// <param name="dbPW">Input argument #3</param>
    /// <returns>An MWArray containing the first output argument.</returns>
    ///
    public MWArray demandForecast(MWArray dbName, MWArray dbID, MWArray dbPW)
    {
      return mcr.EvaluateFunction("demandForecast", dbName, dbID, dbPW);
    }


    /// <summary>
    /// Provides a single output, 4-input MWArrayinterface to the demandForecast MATLAB
    /// function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// parseObj.addParameter('numForecasts',3);
    /// </remarks>
    /// <param name="dbName">Input argument #1</param>
    /// <param name="dbID">Input argument #2</param>
    /// <param name="dbPW">Input argument #3</param>
    /// <param name="VER_ID">Input argument #4</param>
    /// <returns>An MWArray containing the first output argument.</returns>
    ///
    public MWArray demandForecast(MWArray dbName, MWArray dbID, MWArray dbPW, MWArray 
                            VER_ID)
    {
      return mcr.EvaluateFunction("demandForecast", dbName, dbID, dbPW, VER_ID);
    }


    /// <summary>
    /// Provides a single output, 5-input MWArrayinterface to the demandForecast MATLAB
    /// function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// parseObj.addParameter('numForecasts',3);
    /// </remarks>
    /// <param name="dbName">Input argument #1</param>
    /// <param name="dbID">Input argument #2</param>
    /// <param name="dbPW">Input argument #3</param>
    /// <param name="VER_ID">Input argument #4</param>
    /// <param name="periodStart">Input argument #5</param>
    /// <returns>An MWArray containing the first output argument.</returns>
    ///
    public MWArray demandForecast(MWArray dbName, MWArray dbID, MWArray dbPW, MWArray 
                            VER_ID, MWArray periodStart)
    {
      return mcr.EvaluateFunction("demandForecast", dbName, dbID, dbPW, VER_ID, periodStart);
    }


    /// <summary>
    /// Provides a single output, 6-input MWArrayinterface to the demandForecast MATLAB
    /// function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// parseObj.addParameter('numForecasts',3);
    /// </remarks>
    /// <param name="dbName">Input argument #1</param>
    /// <param name="dbID">Input argument #2</param>
    /// <param name="dbPW">Input argument #3</param>
    /// <param name="VER_ID">Input argument #4</param>
    /// <param name="periodStart">Input argument #5</param>
    /// <param name="periodEnd">Input argument #6</param>
    /// <returns>An MWArray containing the first output argument.</returns>
    ///
    public MWArray demandForecast(MWArray dbName, MWArray dbID, MWArray dbPW, MWArray 
                            VER_ID, MWArray periodStart, MWArray periodEnd)
    {
      return mcr.EvaluateFunction("demandForecast", dbName, dbID, dbPW, VER_ID, periodStart, periodEnd);
    }


    /// <summary>
    /// Provides a single output, 7-input MWArrayinterface to the demandForecast MATLAB
    /// function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// parseObj.addParameter('numForecasts',3);
    /// </remarks>
    /// <param name="dbName">Input argument #1</param>
    /// <param name="dbID">Input argument #2</param>
    /// <param name="dbPW">Input argument #3</param>
    /// <param name="VER_ID">Input argument #4</param>
    /// <param name="periodStart">Input argument #5</param>
    /// <param name="periodEnd">Input argument #6</param>
    /// <param name="fcstStart">Input argument #7</param>
    /// <returns>An MWArray containing the first output argument.</returns>
    ///
    public MWArray demandForecast(MWArray dbName, MWArray dbID, MWArray dbPW, MWArray 
                            VER_ID, MWArray periodStart, MWArray periodEnd, MWArray 
                            fcstStart)
    {
      return mcr.EvaluateFunction("demandForecast", dbName, dbID, dbPW, VER_ID, periodStart, periodEnd, fcstStart);
    }


    /// <summary>
    /// Provides a single output, 8-input MWArrayinterface to the demandForecast MATLAB
    /// function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// parseObj.addParameter('numForecasts',3);
    /// </remarks>
    /// <param name="dbName">Input argument #1</param>
    /// <param name="dbID">Input argument #2</param>
    /// <param name="dbPW">Input argument #3</param>
    /// <param name="VER_ID">Input argument #4</param>
    /// <param name="periodStart">Input argument #5</param>
    /// <param name="periodEnd">Input argument #6</param>
    /// <param name="fcstStart">Input argument #7</param>
    /// <param name="fcstEnd">Input argument #8</param>
    /// <returns>An MWArray containing the first output argument.</returns>
    ///
    public MWArray demandForecast(MWArray dbName, MWArray dbID, MWArray dbPW, MWArray 
                            VER_ID, MWArray periodStart, MWArray periodEnd, MWArray 
                            fcstStart, MWArray fcstEnd)
    {
      return mcr.EvaluateFunction("demandForecast", dbName, dbID, dbPW, VER_ID, periodStart, periodEnd, fcstStart, fcstEnd);
    }


    /// <summary>
    /// Provides a single output, 9-input MWArrayinterface to the demandForecast MATLAB
    /// function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// parseObj.addParameter('numForecasts',3);
    /// </remarks>
    /// <param name="dbName">Input argument #1</param>
    /// <param name="dbID">Input argument #2</param>
    /// <param name="dbPW">Input argument #3</param>
    /// <param name="VER_ID">Input argument #4</param>
    /// <param name="periodStart">Input argument #5</param>
    /// <param name="periodEnd">Input argument #6</param>
    /// <param name="fcstStart">Input argument #7</param>
    /// <param name="fcstEnd">Input argument #8</param>
    /// <param name="varargin">Array of MWArrays representing the input arguments 9
    /// through varargin.length+8</param>
    /// <returns>An MWArray containing the first output argument.</returns>
    ///
    public MWArray demandForecast(MWArray dbName, MWArray dbID, MWArray dbPW, MWArray 
                            VER_ID, MWArray periodStart, MWArray periodEnd, MWArray 
                            fcstStart, MWArray fcstEnd, params MWArray[] varargin)
    {
      MWArray[] args = {dbName, dbID, dbPW, VER_ID, periodStart, periodEnd, fcstStart, fcstEnd};
      int nonVarargInputNum = args.Length;
      int varargInputNum = varargin.Length;
      int totalNumArgs = varargInputNum > 0 ? (nonVarargInputNum + varargInputNum) : nonVarargInputNum;
      Array newArr = Array.CreateInstance(typeof(MWArray), totalNumArgs);

      int idx = 0;

      for (idx = 0; idx < nonVarargInputNum; idx++)
        newArr.SetValue(args[idx],idx);

      if (varargInputNum > 0)
      {
        for (int i = 0; i < varargInputNum; i++)
        {
          newArr.SetValue(varargin[i], idx);
          idx++;
        }
      }

      return mcr.EvaluateFunction("demandForecast", (MWArray[])newArr );
    }


    /// <summary>
    /// Provides the standard 0-input MWArray interface to the demandForecast MATLAB
    /// function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// parseObj.addParameter('numForecasts',3);
    /// </remarks>
    /// <param name="numArgsOut">The number of output arguments to return.</param>
    /// <returns>An Array of length "numArgsOut" containing the output
    /// arguments.</returns>
    ///
    public MWArray[] demandForecast(int numArgsOut)
    {
      return mcr.EvaluateFunction(numArgsOut, "demandForecast", new MWArray[]{});
    }


    /// <summary>
    /// Provides the standard 1-input MWArray interface to the demandForecast MATLAB
    /// function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// parseObj.addParameter('numForecasts',3);
    /// </remarks>
    /// <param name="numArgsOut">The number of output arguments to return.</param>
    /// <param name="dbName">Input argument #1</param>
    /// <returns>An Array of length "numArgsOut" containing the output
    /// arguments.</returns>
    ///
    public MWArray[] demandForecast(int numArgsOut, MWArray dbName)
    {
      return mcr.EvaluateFunction(numArgsOut, "demandForecast", dbName);
    }


    /// <summary>
    /// Provides the standard 2-input MWArray interface to the demandForecast MATLAB
    /// function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// parseObj.addParameter('numForecasts',3);
    /// </remarks>
    /// <param name="numArgsOut">The number of output arguments to return.</param>
    /// <param name="dbName">Input argument #1</param>
    /// <param name="dbID">Input argument #2</param>
    /// <returns>An Array of length "numArgsOut" containing the output
    /// arguments.</returns>
    ///
    public MWArray[] demandForecast(int numArgsOut, MWArray dbName, MWArray dbID)
    {
      return mcr.EvaluateFunction(numArgsOut, "demandForecast", dbName, dbID);
    }


    /// <summary>
    /// Provides the standard 3-input MWArray interface to the demandForecast MATLAB
    /// function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// parseObj.addParameter('numForecasts',3);
    /// </remarks>
    /// <param name="numArgsOut">The number of output arguments to return.</param>
    /// <param name="dbName">Input argument #1</param>
    /// <param name="dbID">Input argument #2</param>
    /// <param name="dbPW">Input argument #3</param>
    /// <returns>An Array of length "numArgsOut" containing the output
    /// arguments.</returns>
    ///
    public MWArray[] demandForecast(int numArgsOut, MWArray dbName, MWArray dbID, MWArray 
                              dbPW)
    {
      return mcr.EvaluateFunction(numArgsOut, "demandForecast", dbName, dbID, dbPW);
    }


    /// <summary>
    /// Provides the standard 4-input MWArray interface to the demandForecast MATLAB
    /// function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// parseObj.addParameter('numForecasts',3);
    /// </remarks>
    /// <param name="numArgsOut">The number of output arguments to return.</param>
    /// <param name="dbName">Input argument #1</param>
    /// <param name="dbID">Input argument #2</param>
    /// <param name="dbPW">Input argument #3</param>
    /// <param name="VER_ID">Input argument #4</param>
    /// <returns>An Array of length "numArgsOut" containing the output
    /// arguments.</returns>
    ///
    public MWArray[] demandForecast(int numArgsOut, MWArray dbName, MWArray dbID, MWArray 
                              dbPW, MWArray VER_ID)
    {
      return mcr.EvaluateFunction(numArgsOut, "demandForecast", dbName, dbID, dbPW, VER_ID);
    }


    /// <summary>
    /// Provides the standard 5-input MWArray interface to the demandForecast MATLAB
    /// function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// parseObj.addParameter('numForecasts',3);
    /// </remarks>
    /// <param name="numArgsOut">The number of output arguments to return.</param>
    /// <param name="dbName">Input argument #1</param>
    /// <param name="dbID">Input argument #2</param>
    /// <param name="dbPW">Input argument #3</param>
    /// <param name="VER_ID">Input argument #4</param>
    /// <param name="periodStart">Input argument #5</param>
    /// <returns>An Array of length "numArgsOut" containing the output
    /// arguments.</returns>
    ///
    public MWArray[] demandForecast(int numArgsOut, MWArray dbName, MWArray dbID, MWArray 
                              dbPW, MWArray VER_ID, MWArray periodStart)
    {
      return mcr.EvaluateFunction(numArgsOut, "demandForecast", dbName, dbID, dbPW, VER_ID, periodStart);
    }


    /// <summary>
    /// Provides the standard 6-input MWArray interface to the demandForecast MATLAB
    /// function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// parseObj.addParameter('numForecasts',3);
    /// </remarks>
    /// <param name="numArgsOut">The number of output arguments to return.</param>
    /// <param name="dbName">Input argument #1</param>
    /// <param name="dbID">Input argument #2</param>
    /// <param name="dbPW">Input argument #3</param>
    /// <param name="VER_ID">Input argument #4</param>
    /// <param name="periodStart">Input argument #5</param>
    /// <param name="periodEnd">Input argument #6</param>
    /// <returns>An Array of length "numArgsOut" containing the output
    /// arguments.</returns>
    ///
    public MWArray[] demandForecast(int numArgsOut, MWArray dbName, MWArray dbID, MWArray 
                              dbPW, MWArray VER_ID, MWArray periodStart, MWArray 
                              periodEnd)
    {
      return mcr.EvaluateFunction(numArgsOut, "demandForecast", dbName, dbID, dbPW, VER_ID, periodStart, periodEnd);
    }


    /// <summary>
    /// Provides the standard 7-input MWArray interface to the demandForecast MATLAB
    /// function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// parseObj.addParameter('numForecasts',3);
    /// </remarks>
    /// <param name="numArgsOut">The number of output arguments to return.</param>
    /// <param name="dbName">Input argument #1</param>
    /// <param name="dbID">Input argument #2</param>
    /// <param name="dbPW">Input argument #3</param>
    /// <param name="VER_ID">Input argument #4</param>
    /// <param name="periodStart">Input argument #5</param>
    /// <param name="periodEnd">Input argument #6</param>
    /// <param name="fcstStart">Input argument #7</param>
    /// <returns>An Array of length "numArgsOut" containing the output
    /// arguments.</returns>
    ///
    public MWArray[] demandForecast(int numArgsOut, MWArray dbName, MWArray dbID, MWArray 
                              dbPW, MWArray VER_ID, MWArray periodStart, MWArray 
                              periodEnd, MWArray fcstStart)
    {
      return mcr.EvaluateFunction(numArgsOut, "demandForecast", dbName, dbID, dbPW, VER_ID, periodStart, periodEnd, fcstStart);
    }


    /// <summary>
    /// Provides the standard 8-input MWArray interface to the demandForecast MATLAB
    /// function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// parseObj.addParameter('numForecasts',3);
    /// </remarks>
    /// <param name="numArgsOut">The number of output arguments to return.</param>
    /// <param name="dbName">Input argument #1</param>
    /// <param name="dbID">Input argument #2</param>
    /// <param name="dbPW">Input argument #3</param>
    /// <param name="VER_ID">Input argument #4</param>
    /// <param name="periodStart">Input argument #5</param>
    /// <param name="periodEnd">Input argument #6</param>
    /// <param name="fcstStart">Input argument #7</param>
    /// <param name="fcstEnd">Input argument #8</param>
    /// <returns>An Array of length "numArgsOut" containing the output
    /// arguments.</returns>
    ///
    public MWArray[] demandForecast(int numArgsOut, MWArray dbName, MWArray dbID, MWArray 
                              dbPW, MWArray VER_ID, MWArray periodStart, MWArray 
                              periodEnd, MWArray fcstStart, MWArray fcstEnd)
    {
      return mcr.EvaluateFunction(numArgsOut, "demandForecast", dbName, dbID, dbPW, VER_ID, periodStart, periodEnd, fcstStart, fcstEnd);
    }


    /// <summary>
    /// Provides the standard 9-input MWArray interface to the demandForecast MATLAB
    /// function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// parseObj.addParameter('numForecasts',3);
    /// </remarks>
    /// <param name="numArgsOut">The number of output arguments to return.</param>
    /// <param name="dbName">Input argument #1</param>
    /// <param name="dbID">Input argument #2</param>
    /// <param name="dbPW">Input argument #3</param>
    /// <param name="VER_ID">Input argument #4</param>
    /// <param name="periodStart">Input argument #5</param>
    /// <param name="periodEnd">Input argument #6</param>
    /// <param name="fcstStart">Input argument #7</param>
    /// <param name="fcstEnd">Input argument #8</param>
    /// <param name="varargin">Array of MWArrays representing the input arguments 9
    /// through varargin.length+8</param>
    /// <returns>An Array of length "numArgsOut" containing the output
    /// arguments.</returns>
    ///
    public MWArray[] demandForecast(int numArgsOut, MWArray dbName, MWArray dbID, MWArray 
                              dbPW, MWArray VER_ID, MWArray periodStart, MWArray 
                              periodEnd, MWArray fcstStart, MWArray fcstEnd, params 
                              MWArray[] varargin)
    {
      MWArray[] args = {dbName, dbID, dbPW, VER_ID, periodStart, periodEnd, fcstStart, fcstEnd};
      int nonVarargInputNum = args.Length;
      int varargInputNum = varargin.Length;
      int totalNumArgs = varargInputNum > 0 ? (nonVarargInputNum + varargInputNum) : nonVarargInputNum;
      Array newArr = Array.CreateInstance(typeof(MWArray), totalNumArgs);

      int idx = 0;

      for (idx = 0; idx < nonVarargInputNum; idx++)
        newArr.SetValue(args[idx],idx);

      if (varargInputNum > 0)
      {
        for (int i = 0; i < varargInputNum; i++)
        {
          newArr.SetValue(varargin[i], idx);
          idx++;
        }
      }

      return mcr.EvaluateFunction(numArgsOut, "demandForecast", (MWArray[])newArr );
    }


    /// <summary>
    /// Provides an interface for the demandForecast function in which the input and
    /// output
    /// arguments are specified as an array of MWArrays.
    /// </summary>
    /// <remarks>
    /// This method will allocate and return by reference the output argument
    /// array.<newpara></newpara>
    /// M-Documentation:
    /// parseObj.addParameter('numForecasts',3);
    /// </remarks>
    /// <param name="numArgsOut">The number of output arguments to return</param>
    /// <param name= "argsOut">Array of MWArray output arguments</param>
    /// <param name= "argsIn">Array of MWArray input arguments</param>
    ///
    public void demandForecast(int numArgsOut, ref MWArray[] argsOut, MWArray[] argsIn)
    {
      mcr.EvaluateFunction("demandForecast", numArgsOut, ref argsOut, argsIn);
    }



    /// <summary>
    /// This method will cause a MATLAB figure window to behave as a modal dialog box.
    /// The method will not return until all the figure windows associated with this
    /// component have been closed.
    /// </summary>
    /// <remarks>
    /// An application should only call this method when required to keep the
    /// MATLAB figure window from disappearing.  Other techniques, such as calling
    /// Console.ReadLine() from the application should be considered where
    /// possible.</remarks>
    ///
    public void WaitForFiguresToDie()
    {
      mcr.WaitForFiguresToDie();
    }



    #endregion Methods

    #region Class Members

    private static MWMCR mcr= null;

    private static Exception ex_= null;

    private bool disposed= false;

    #endregion Class Members
  }
}
