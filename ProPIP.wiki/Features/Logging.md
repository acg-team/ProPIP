[Home](../ProPIP/-Progressive-Multiple-Sequence-Alignment-with-Poisson-Indel-Process.md) | [Features](Index.md) |

---

### Log files verbosity

ProPIP takes advantage of the Google Log libraries to generate debugging informations about the program execution. To tune the options according to the need of the user here a short list of the options:

#### 1. Output device

    export GLOG_logtostderr (bool, default=false)
        Log messages to stderr instead of logfiles.
        Note: you can set binary flags to true by specifying 1, true, or yes (case insensitive). Also, you can set binary flags to false by specifying 0, false, or no (again, case insensitive).

#### 2. Threshold

    export GLOG_stderrthreshold (int, default=2, which is ERROR)
        Copy log messages at or above this level to stderr in addition to logfiles. The numbers of severity levels INFO, WARNING, ERROR, and FATAL are 0, 1, 2, and 3, respectively.

#### 3. Minimum level

    export GLOG_minloglevel (int, default=0, which is INFO)
        Log messages at or above this level. Again, the numbers of severity levels INFO, WARNING, ERROR, and FATAL are 0, 1, 2, and 3, respectively.

#### 4. Location of the log files

    export GLOG_log_dir (string, default="")
        If specified, logfiles are written into this directory instead of the default logging directory.

#### 5. Level of the verbose output

    export GLOG_v (int, default=0, min=0, max=3)
        Show all VLOG(m) messages for m less or equal the value of this flag. Overridable by --vmodule. See the section about verbose logging for more detail.

#### 6. Module verbosity

    export GLOG_vmodule (string, default="")
        Per-module verbose level. The argument has to contain a comma-separated list of <module name>=<log level>. <module name> is a glob pattern (e.g., gfs* for all modules whose name starts with "gfs"), matched against the filename base (that is, name ignoring .cc/.h./-inl.h). <log level> overrides any value given by --v. See also the section about verbose logging.



### Example

Module verbosity levels:

    export GLOG_vmodule=TreeRearrangment=2,UnifiedTSHTopologySearch=2,TreeSearch=2

Redirect LOG to stderr:
    
    export GLOG_logtostderr=1
